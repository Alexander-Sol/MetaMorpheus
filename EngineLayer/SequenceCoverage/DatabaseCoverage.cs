using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.IO;
using UsefulProteomicsDatabases;

namespace EngineLayer.SequenceCoverage
{
    internal class DatabaseCoverage
    {
        public List<Protein> ProteinList { get; private set; }
        public List<PsmFromTsv> ReadPsms { get; private set; }
        public Dictionary<Protein, ProteinCoverage> CoverageDictionary { get; private set; }

        public DatabaseCoverage(string databasePath, string peptidePath)
        {
            ProteinList = LoadProteins(databasePath);
            ReadPsms = PsmTsvReader.ReadTsv(peptidePath, out var warnings);
            foreach (var protein in ProteinList)
            {
                var proteinPeptides = ReadPsms.Where(p => p.ProteinAccession.Equals(protein.Accession)).ToList();
                CoverageDictionary.Add(protein, new ProteinCoverage(protein, proteinPeptides));
            }
        }
        protected List<Protein> LoadProteins(string databasePath)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(databasePath).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(databasePath)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                proteinList = ProteinDbLoader.LoadProteinFasta(databasePath, true, DecoyType.None, false, out dbErrors, ProteinDbLoader.UniprotAccessionRegex,
                    ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);
                if (!proteinList.Any())
                {
                    //Warn("Warning: No protein entries were found in the database");
                    return new List<Protein>() { };
                }
                else
                {
                    return proteinList;
                }

            }
            else
            {
                List<string> modTypesToExclude = new List<string> { };
                proteinList = ProteinDbLoader.LoadProteinXML(databasePath, true, DecoyType.None, GlobalVariables.AllModsKnown, false, modTypesToExclude,
                    out Dictionary<string, Modification> um, -1, 4, 1);
                if (!proteinList.Any())
                {
                    //Warn("Warning: No protein entries were found in the database");
                    return new List<Protein>() { };
                }
                else
                {
                    return proteinList;
                }
            }
        }
    }
}

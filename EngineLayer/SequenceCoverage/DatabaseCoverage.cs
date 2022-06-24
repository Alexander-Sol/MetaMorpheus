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
    public class DatabaseCoverage
    {
        public List<Protein> ProteinList { get; private set; }
        public List<PsmFromTsv> ReadPsms { get; private set; }
        private Dictionary<Protein, ProteinCoverage> CoverageDictionary { get; set; }
        public List<string> ProteasesPresent { get; set; }
        public string OutputFolder { get; private set; }

        public DatabaseCoverage(string databasePath, string peptidePath)
        {
            ProteinList = LoadProteins(databasePath);
            ReadPsms = PsmTsvReader.ReadTsv(peptidePath, out var warnings);
            IEnumerable<string> previousAAs = ReadPsms.Select(p => p.PreviousAminoAcid).Distinct();
            ProteasesPresent = new List<string>();
            if (previousAAs.Contains("D") | previousAAs.Contains("E")) ProteasesPresent.Add("gluC");
            if (previousAAs.Contains("P") | previousAAs.Contains("A")) ProteasesPresent.Add("proA");
            if (previousAAs.Contains("K") | previousAAs.Contains("R")) ProteasesPresent.Add("trypsin");

            CoverageDictionary = new();
            foreach (var protein in ProteinList)
            {
                var proteinPeptides = ReadPsms
                    .Where(p => p.PEP_QValue < 0.1)
                    .Where(p => p.ProteinAccession.Equals(protein.Accession))
                    .ToList();
                CoverageDictionary.Add(protein, new ProteinCoverage(protein, proteinPeptides));
            }
        }
        public void WriteCoverageResults(string outputFolder)
        {
            OutputFolder = outputFolder;
            string outputPath = Path.Combine(outputFolder, @"CoverageResults.tsv");
            using (var output = new StreamWriter(outputPath))
            {
                output.WriteLine(DatabaseCoverage.GetTabSeparatedHeader());
                foreach (ProteinCoverage protein in CoverageDictionary.Values)
                {
                    output.WriteLine(protein.ToString());
                }
            }
        }
        public static string GetTabSeparatedHeader()
        {
            List<string> sb = new();
            sb.Add("Protein Accesion");
            sb.Add("Total Peptide Coverage");
            sb.Add("Total Fragment Coverage");
            sb.Add("gluC Peptide Coverage");
            sb.Add("gluC Fragment Coverage");
            sb.Add("proA Peptide Coverage");
            sb.Add("proA Fragment Coverage");
            sb.Add("trypsin Peptide Coverage");
            sb.Add("trypsin Fragment Coverage");
            string header = String.Join('\t', sb);
            return header;
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

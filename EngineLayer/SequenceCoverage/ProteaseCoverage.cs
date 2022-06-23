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
    internal class ProteaseCoverage
    {
        public readonly Protease Protease;
        public readonly Protein Protein;
        public readonly Char[] Residues;
        public bool[] PeptideCoverage { get; private set; }
        public int[] FragmentCoverage { get; private set; }

        public ProteaseCoverage()
        {

        }

        protected List<Protein> LoadProteins(DbForDigestion database)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(database.FilePath).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(database.FilePath)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                proteinList = ProteinDbLoader.LoadProteinFasta(database.FilePath, true, DecoyType.None, false, out dbErrors, ProteinDbLoader.UniprotAccessionRegex,
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
                proteinList = ProteinDbLoader.LoadProteinXML(database.FilePath, true, DecoyType.None, GlobalVariables.AllModsKnown, false, modTypesToExclude,
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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Proteomics.Fragmentation;
using System.IO;
using UsefulProteomicsDatabases;

namespace EngineLayer.SequenceCoverage
{
    internal class ProteinCoverage
    {
        public readonly Char[] Residues;
        public readonly Protein Protein;
        public bool[] PeptideCoverage { get; private set; }
        public bool[] FragmentCoverage { get; private set; }
        public int[] PepBondCoverage { get; private set; }
        public Dictionary<string, ProteaseCoverage> ProteaseCoverageDict { get; private set; } //At some point, need to incorporate actual protease as the key for this dict
        public float PepCoveragePercent { get; private set; }
        public float FragmentCoveragePercent { get; private set; }

        public static Dictionary<string, List<string>> proteaseAminoAcidPairs = new Dictionary<string, List<string>>
        {
            { "gluC", new List<string> { "E", "D", "M", "-"} },
            { "proA", new List<string> { "P", "A", "M", "-"} },
            { "trypsin", new List<string> { "K", "R", "M", "-"} },
        };

        public ProteinCoverage(Protein protein, List<PsmFromTsv> readPsms )
        {
            Residues = protein.BaseSequence.ToCharArray();
            Protein = protein;
            PeptideCoverage = new bool[Protein.Length];     // An array representing every residue within a given protein
            PepBondCoverage = new int[Protein.Length - 1]; // An array representing all the peptide bonds within a given protein
            FragmentCoverage = new bool[Protein.Length];
            ProteaseCoverageDict = new();

            // This is a shitty manual work around. At some point, need to incorporate actual protease objects
            // However, psmtsv don't store the protease used to generate the peptide. 
            List<PsmFromTsv> gluCPsms = readPsms
                .Where(p => proteaseAminoAcidPairs["gluC"].Contains(p.PreviousAminoAcid))
                .ToList();
            List<PsmFromTsv> proAPsms = readPsms
                .Where(p => proteaseAminoAcidPairs["proA"].Contains(p.PreviousAminoAcid))
                .ToList();
            List<PsmFromTsv> trypsinPsms = readPsms
                .Where(p => proteaseAminoAcidPairs["trypsin"].Contains(p.PreviousAminoAcid))
                .ToList();

            ProteaseCoverageDict.Add("gluC", new ProteaseCoverage("gluC", Protein, Residues, gluCPsms));
            ProteaseCoverageDict.Add("proA", new ProteaseCoverage("proA", Protein, Residues, proAPsms));
            ProteaseCoverageDict.Add("trypsin", new ProteaseCoverage("trypsin", Protein, Residues, trypsinPsms));

            DetermineTotalPeptideCoverage(ProteaseCoverageDict, PeptideCoverage);
            DetermineTotalPepBondCoverage(ProteaseCoverageDict, PepBondCoverage);
            DetermineFragmentCoverage(PepBondCoverage, FragmentCoverage);
            CalculatePercentCoverage();
            var place = 0;
        }

        public void DetermineTotalPeptideCoverage(Dictionary<string, ProteaseCoverage> proteaseCoverageDict, bool[] peptideCoverage)
        {
            bool[][] combinedCoverage = new bool[proteaseCoverageDict.Count][];
            int column = 0;
            foreach(ProteaseCoverage proteaseCov in proteaseCoverageDict.Values)
            {
                combinedCoverage[column] = proteaseCov.PeptideCoverage;
                column++;
            }
            for (int i = 0; i < peptideCoverage.Length; i++)
            {
                for(int j = 0; j < proteaseCoverageDict.Count; j++)
                {
                    peptideCoverage[i] = peptideCoverage[i] | combinedCoverage[j][i];
                }
            }
        }

        public void DetermineTotalPepBondCoverage(Dictionary<string, ProteaseCoverage> proteaseCoverageDict, int[] pepBondCoverage)
        {
            int[][] combinedCoverage = new int[proteaseCoverageDict.Count][];
            int column = 0;
            foreach (ProteaseCoverage proteaseCov in proteaseCoverageDict.Values)
            {
                combinedCoverage[column] = proteaseCov.PepBondCoverage;
                column++;
            }
            for (int i = 0; i < pepBondCoverage.Length; i++)
            {
                for (int j = 0; j < proteaseCoverageDict.Count; j++)
                {
                    pepBondCoverage[i] += combinedCoverage[j][i];
                }
            }
        }

        public void DetermineFragmentCoverage(int[] pepBondCoverage, bool[] fragmentCoverage)
        {
            int n;
            int nplus1;
            // If sequential peptide bonds were fragmented, then the enclosed residue is marked as being covered (within fragmentCoverage)
            for(int i = -1; i < pepBondCoverage.Length; i++)
            {
                n = i < 0 ? 1 : pepBondCoverage[i];
                nplus1 = (i+1) == pepBondCoverage.Length ? 1 : pepBondCoverage[i + 1];
                fragmentCoverage[i + 1] = n * nplus1 > 0 ? true : false;
            }
        }

        public void CalculatePercentCoverage()
        {
            PepCoveragePercent = (float)100.00 * PeptideCoverage.Count(b => b) / Protein.Length;
            FragmentCoveragePercent = (float)100.00 * FragmentCoverage.Count(b => b) / Protein.Length;
        }

        public (float, float) CombinatorialCoverage(ProteaseCoverage a, ProteaseCoverage b)
        {
            return ( (float)0.0, (float)0.0 );
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, string>());
        }

        private string ToString(Dictionary<string, string> CoverageWriter)
        {
            CoverageWriter.Add("Protein Accession", Protein.Accession);
            CoverageWriter.Add("Total Peptide Coverage", PepCoveragePercent.ToString("0.00"));
            CoverageWriter.Add("Total Fragment Coverage", FragmentCoveragePercent.ToString("0.00"));
            foreach(ProteaseCoverage proteaseCoverage in ProteaseCoverageDict.Values)
            {
                CoverageWriter.Add(proteaseCoverage.Protease + " Peptide Coverage", proteaseCoverage.PepCoveragePercent.ToString("0.00"));
                CoverageWriter.Add(proteaseCoverage.Protease + " Fragment Coverage", proteaseCoverage.FragmentCoveragePercent.ToString("0.00"));
            }
            return string.Join("\t", CoverageWriter.Values);
        }
    }
}

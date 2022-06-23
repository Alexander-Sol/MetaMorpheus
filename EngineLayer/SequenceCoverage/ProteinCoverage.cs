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
        public ProteinCoverage(Protein protein, List<PsmFromTsv> readPsms )
        {
            Residues = protein.BaseSequence.ToCharArray();
            Protein = protein;
            PeptideCoverage = new bool[Protein.Length];     // An array representing every residue within a given protein
            FragmentCoverage = new bool[Protein.Length];
            PepBondCoverage = new int[Protein.Length - 1]; // An array representing all the peptide bonds within a given protein

            DeterminePeptideCoverage(PeptideCoverage, readPsms);
            DeterminePepBondCoverage(PepBondCoverage, readPsms);
            DetermineFragmentCoverage(PepBondCoverage, FragmentCoverage);

            var place = 0;
        }

        public void DeterminePeptideCoverage(bool[] coverageArray, List<PsmFromTsv> readPsms)
        {
            int startResidue; 
            int endResidue;
            foreach(PsmFromTsv peptide in readPsms)
            {
                startResidue = Int32.Parse(peptide.StartResidueInProtein);
                endResidue = Int32.Parse(peptide.EndResidueInProtein);
                for (int i = startResidue - 1; i < endResidue; i++) coverageArray[i] = true;
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

        public void DeterminePepBondCoverage(int[] coverageArray, List<PsmFromTsv> readPsms)
        {
            int startResidue;
            int endResidue;
            foreach (PsmFromTsv peptide in readPsms)
            {
                startResidue = Int32.Parse(peptide.StartResidueInProtein);
                endResidue = Int32.Parse(peptide.EndResidueInProtein);
                //N-terminal
                List<MatchedFragmentIon> nTermFragments = peptide.MatchedIons.Where(x => x.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N).ToList();
                //C-terminal in reverse order
                List<MatchedFragmentIon> cTermFragments = peptide.MatchedIons.Where(x => x.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C).ToList();

                foreach (MatchedFragmentIon ion in nTermFragments)
                {
                    coverageArray[startResidue + ion.NeutralTheoreticalProduct.AminoAcidPosition - 2]++;
                }
                foreach (MatchedFragmentIon ion in cTermFragments)
                {
                    coverageArray[startResidue + ion.NeutralTheoreticalProduct.AminoAcidPosition - 2]++;
                }
            }
        }
    }
}

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
        public int[] FragmentCoverage { get; private set; }
        public ProteinCoverage(Protein protein, List<PsmFromTsv> readPsms )
        {
            Residues = protein.BaseSequence.ToCharArray();
            Protein = protein;
            PeptideCoverage = new bool[Protein.Length];
            FragmentCoverage = new int[Protein.Length];

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

        public void DetermineFragmentCoverage(int[] coverageArray, List<PsmFromTsv> readPsms)
        {
            foreach (PsmFromTsv peptide in readPsms)
            {
                //N-terminal
                List<MatchedFragmentIon> nTermFragments = peptide.MatchedIons.Where(x => x.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N).ToList();
                //C-terminal in reverse order
                List<MatchedFragmentIon> cTermFragments = peptide.MatchedIons.Where(x => x.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C).OrderByDescending(x => x.NeutralTheoreticalProduct.FragmentNumber).ToList();

                int position = nTermFragments[0].NeutralTheoreticalProduct.AminoAcidPosition;
                int place = 0;
            }
        }
    }
}

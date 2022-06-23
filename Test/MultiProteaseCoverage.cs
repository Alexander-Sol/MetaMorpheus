using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class MultiProteaseCoverage
    {
        [Test]
        public static void MultiProteaseTest()
        {
            string psmFile = @"C:\Users\Alex\Desktop\MutantProteinSP2022\EnterolysinA\ProA_EThcD_Search\Task1-SearchTask\AllPeptides.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Dictionary<string, List<PsmFromTsv>> peptidesByProtein = new Dictionary<string, List<PsmFromTsv>>();
            string[] uniqueProteins = parsedPsms.Select(p => p.ProteinAccession).Distinct().ToArray();
            foreach (string protein in uniqueProteins)
            {
                List<PsmFromTsv> psmForProtein = parsedPsms.Where(p => p.ProteinAccession.Equals(protein)).ToList();
                peptidesByProtein.Add(protein, psmForProtein);
            }

            Dictionary<string, bool[]> peptideLevelCoverage = new Dictionary<string, bool[]>();
            foreach (string protein in peptidesByProtein.Keys)
            {
                int numberOfResidues = peptidesByProtein[protein].Select(p => Int32.Parse(p.EndResidueInProtein)).Max();
                bool[] pepFound = new bool[numberOfResidues]; //creates array, one bool for each residue. Initially, every entry is false
                int startResidue;
                int endResidue;
                foreach (PsmFromTsv peptide in peptidesByProtein[protein])
                {
                    startResidue = Int32.Parse(peptide.StartResidueInProtein);
                    endResidue = Int32.Parse(peptide.EndResidueInProtein);
                    for(int i = startResidue - 1; i < endResidue; i++) pepFound[i] = true;
                }

                peptideLevelCoverage.Add(protein, pepFound);
                Console.WriteLine(protein + " has " + 100 * pepFound.Count(c => c) / numberOfResidues + "% sequence coverage");
            }

            

            Assert.AreEqual(2, parsedPsms.Count);
        }

    }
}

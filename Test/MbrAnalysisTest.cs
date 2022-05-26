﻿using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using EngineLayer;
using FlashLFQ;
using MathNet.Numerics.Statistics; // Necessary for calculating correlation 
using System.Text.RegularExpressions;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using EngineLayer.ClassicSearch;


namespace Test
{
    [TestFixture]
    public static class
    MbrAnalysisTest
    {
        [Test]
        public static void MbrPostSearchAnalysisTest()
        {
            // This block of code converts from PsmFromTsv to PeptdieSpectralMatch objects.
            // It also deals with one specific Carbamidomethylation, defined in advance
            string psmtsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MsMsids.psmtsv");
            List<PsmFromTsv> tsvPsms = PsmTsvReader.ReadTsv(psmtsvPath, out var warnings);
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();
            List<Protein> proteinList = new List<Protein>();
            MyFileManager myFileManager = new MyFileManager(true);

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod1 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            Dictionary<string, Modification> carbamidoDict = new Dictionary<string, Modification> { { "Carbamidomethyl on C", mod1 } };
            
            foreach (PsmFromTsv readPsm in tsvPsms)
            {
                string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "TestData", "MbrAnalysisTest", readPsm.FileNameWithoutExtension + ".mzML");
                MsDataScan scan = myFileManager.LoadFile(filePath, new CommonParameters()).GetOneBasedScan(readPsm.Ms2ScanNumber);
                Ms2ScanWithSpecificMass ms2Scan = new Ms2ScanWithSpecificMass(scan, readPsm.PrecursorMz, readPsm.PrecursorCharge,
                    filePath, new CommonParameters());
                Protein protein = new Protein(readPsm.BaseSeq, readPsm.ProteinAccession, readPsm.OrganismName);
                string[] startAndEndResidues = readPsm.StartAndEndResiduesInProtein.Split(" ");
                int startResidue = Int32.Parse(startAndEndResidues[0].Trim('['));
                int endResidue = Int32.Parse(startAndEndResidues[2].Trim(']'));
                PeptideWithSetModifications pwsm = readPsm.FullSequence.Contains("[") ?
                    new PeptideWithSetModifications(readPsm.FullSequence, carbamidoDict, p: protein, digestionParams: new DigestionParams(), oneBasedStartResidueInProtein: startResidue, oneBasedEndResidueInProtein: endResidue) :
                    new PeptideWithSetModifications(readPsm.FullSequence, null, p: protein, digestionParams: new DigestionParams(), oneBasedStartResidueInProtein: startResidue, oneBasedEndResidueInProtein: endResidue);
                PeptideSpectralMatch psm = new PeptideSpectralMatch(pwsm, 0, readPsm.Score, readPsm.Ms2ScanNumber, ms2Scan,
                    new CommonParameters(), readPsm.MatchedIons);
                
                psms.Add(psm);
                proteinList.Add(protein);
            }

            List<string> rawSlices = new List<string> {
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MbrTest_J3.mzML"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MbrTest_K13.mzML") };
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput/individual"));
            Dictionary<string, int[]> numSpectraPerFile = new Dictionary<string, int[]> { { "MbrTest_J3", new int[] { 8, 8 } }, { "MbrTest_K13", new int[] { 8, 8 } } };
            List<DbForTask> databaseList = new List<DbForTask>() {new DbForTask(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\HumanFastaSlice.fasta"), false) };
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput");

            SearchTask searchTask = new SearchTask
            {
                SearchParameters = new SearchParameters()
                {
                    DoQuantification = false,
                    WriteSpectralLibrary = false,
                    MatchBetweenRuns = false,
                    DoMbrAnalysis = false,
                    WriteMzId = false
                },
                CommonParameters = new CommonParameters()
            };

            var testTaskResults = searchTask.RunTask(outputFolder, databaseList, rawSlices, "name");

            PostSearchAnalysisTask postSearchTask = new PostSearchAnalysisTask()
            {
                Parameters = new PostSearchAnalysisParameters()
                {
                    ProteinList = proteinList,
                    AllPsms = psms,
                    CurrentRawFileList = rawSlices,
                    DatabaseFilenameList = databaseList,
                    OutputFolder = outputFolder,
                    NumMs2SpectraPerFile = numSpectraPerFile,
                    ListOfDigestionParams = new HashSet<DigestionParams> { new DigestionParams(generateUnlabeledProteinsForSilac: false) },
                    SearchTaskResults = testTaskResults,
                    MyFileManager = myFileManager,
                    IndividualResultsOutputFolder = Path.Combine(outputFolder, "individual"),
                    SearchParameters = new SearchParameters()
                    {
                        DoQuantification = true,
                        WriteSpectralLibrary = true,
                        MatchBetweenRuns = true,
                        DoMbrAnalysis = true,
                        WriteMzId = false,
                        QuantifyPpmTol = 25
                    }
                },
                CommonParameters = new CommonParameters(),
                FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> {
                    (rawSlices[0], new CommonParameters()),
                    (rawSlices[1], new CommonParameters()) 
                }
            };

            int placeholder = 1;

            postSearchTask.Run();

            placeholder = 0;

            Assert.That(placeholder == 0);

        }

        [Test]
        public static void MiniClassicSearchEngineTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.AreEqual(7, allPsmsArray.Length);
            Assert.IsTrue(allPsmsArray[5].Score > 38);
            Assert.AreEqual("VIHDNFGIVEGLMTTVHAITATQK", allPsmsArray[5].BaseSequence);
            Assert.IsTrue(!allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));


            foreach (PeptideSpectralMatch psm in allPsmsArray.Where(p => p != null))
            {
                PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                List<(string FileName, CommonParameters Parameters)> fileSpecificParameters = new() { (myFile, commonParameters) };
                PeptideSpectralMatch[] peptideSpectralMatches = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

                new MiniClassicSearchEngine(pwsm, peptideSpectralMatches, listOfSortedms2Scans, variableModifications, fixedModifications, searchModes,
                    commonParameters, fileSpecificParameters, sl, new List<string>()).Run();

                Assert.AreEqual(allPsmsArray[5].BaseSequence, peptideSpectralMatches[5].BaseSequence);
                Assert.That(peptideSpectralMatches[5].SpectralAngle, Is.EqualTo(allPsmsArray[5].SpectralAngle).Within(0.01));
            }
        }
    }
}
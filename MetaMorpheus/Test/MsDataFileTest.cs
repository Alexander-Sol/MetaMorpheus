using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Readers;
using TaskLayer;
using MsDataFile = MassSpectrometry.MsDataFile;

namespace Test
{
    [TestFixture]
    public static class MsDataFileTest
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public static void GetIsolationWindows()
        {
            List<string> jurkatFiles = new List<string>
            {
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_6.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_1.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_2.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_3.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_4.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_1\20100614_Velos1_TaGe_SA_Jurkat_5.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_01.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_02.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_03.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_04.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_05.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_2\20100730_Velos1_TaGe_SA_Jurkat_06_100731121305.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat5.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat6.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat1.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat2.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat3.raw",
                @"D:\JurkatTrypsin\Jurkat_Mann11\Jurkat_3\20101230_Velos1_TaGe_SA_Jurkat4.raw"

            };

            string[] header = new string[]
            {
                "File",
                "FileScan",
                "Scan Number",
                "Isolation Center",
                "Isolation Width"
            };
            using StreamWriter writer = new StreamWriter(@"D:\JurkatTrypsin\Jurkat_Mann11\IsolationTable.tsv");
            writer.WriteLine(String.Join('\t', header));

            foreach (string jurkatFile in jurkatFiles)
            {
                MsDataFile dataFile = MsDataFileReader.GetDataFile(jurkatFile);
                var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2);
                string fileName = Path.GetFileName(jurkatFile);

                foreach (var scan in ms2Scans)
                {
                    string scanNumber = scan.OneBasedScanNumber.ToString();
                    string[] scanInfo = new string[]
                    {
                        fileName,
                        fileName + scanNumber,
                        scanNumber,
                        scan.IsolationMz.ToString(),
                        scan.IsolationWidth.ToString()
                    };
                    writer.WriteLine(String.Join('\t', scanInfo));
                }
            }
        }

        [Test]
        public static void TestLoadAndRunMgf()
        {
            //The purpose of this test is to ensure that mgfs can be run without crashing.
            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string mgfName = @"TestData\ok.mgf";
            string xmlName = @"TestData\okk.xml";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestLoadAndRunMgf");

            SearchTask task1 = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = true
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new()
            {
                ("task1", task1),
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestCompressionDecompression()
        {
            string testInputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CompressionTest");
            DirectoryInfo testDirectory = new(testInputFolder);
            MyFileManager.CompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreEqual(".gz", file.Extension);
            }

            MyFileManager.DecompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreNotEqual(".gz", file.Extension);
            }
        }

        [Test]
        public static void TestMs2ScanWithSpecificMass()
        {
            Ms2ScanWithSpecificMass scanB = new(
                new MsDataScan(
                    new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            var closestExperimentalMassB = scanB.GetClosestExperimentalIsotopicEnvelope(10);

            Assert.IsNull(closestExperimentalMassB);
        }
    }
}
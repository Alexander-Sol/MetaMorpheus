using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class IonQuantReader
    {
        public static readonly Dictionary<string, string> SampleToFileDict =
            new Dictionary<string, string>
            {
                { "nanoPOTS_1", @"D:\HelaSingleCellQCmzML\Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML" },
                { "nanoPOTS_2", @"D:\HelaSingleCellQCmzML\Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML" },
                { "nanoPOTS_3", @"D:\HelaSingleCellQCmzML\Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML" },
                { "nanoPOTS_4", @"D:\HelaSingleCellQCmzML\Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML" },
                { "nanoPOTS_5", @"D:\HelaSingleCellQCmzML\Ex_Auto_K13_30umTB_02ngQC_60m_1.mzML" },
                { "nanoPOTS_6", @"D:\HelaSingleCellQCmzML\Ex_Auto_K13_30umTB_02ngQC_60m_2.mzML" },
                { "nanoPOTS_7", @"D:\HelaSingleCellQCmzML\Ex_Auto_K13_30umTB_2ngQC_60m_1.mzML" },
                { "nanoPOTS_8", @"D:\HelaSingleCellQCmzML\Ex_Auto_K13_30umTB_2ngQC_60m_2.mzML" },
                { "nanoPOTS_9", @"D:\HelaSingleCellQCmzML\Ex_Auto_W17_30umTB_02ngQC_60m_3.mzML" },
                { "nanoPOTS_10", @"D:\HelaSingleCellQCmzML\Ex_Auto_W17_30umTB_02ngQC_60m_4.mzML" },
                { "nanoPOTS_11", @"D:\HelaSingleCellQCmzML\Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML" },
                { "nanoPOTS_12", @"D:\HelaSingleCellQCmzML\Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML" }
            };

        private static MyTaskResults SearchTaskResults;
        private static List<PsmFromTsv> TsvPsms;
        private static List<PeptideSpectralMatch> Psms;
        private static List<Protein> ProteinList;
        private static MyFileManager MyFileManager;
        private static List<string> RawSlices;
        private static List<DbForTask> DatabaseList;
        private static string OutputFolder;
        private static Dictionary<string, int[]> NumSpectraPerFile;
        private Dictionary<string, int> CombinedIonHeader;



        [Test]
        public void ReadCombinedIons()
        {
            StreamReader reader = new(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData",
                @"MbrAnalysisTest\IonQuant\combined_ion.tsv"));
            List<IonQuantInfo> ionQuantMbrPeptides = new();
            CombinedIonHeader = new();

            using (reader)
            {
                string[] header = reader.ReadLine().Split('\t');
                for (int i = 0; i < header.Length; i++)
                {
                    CombinedIonHeader.Add(header[i], i);
                }
    
                while (!reader.EndOfStream)
                {
                    string? line = reader.ReadLine();
                    string[] combinedIonLine = line.Split('\t');
                    foreach (IonQuantInfo ion in ParseCombinedIonsFile(combinedIonLine))
                    {
                        ionQuantMbrPeptides.Add(ion);
                    }
                }

            }

            UpdateIons(ionQuantMbrPeptides);

            // This block of code converts from PsmFromTsv to PeptideSpectralMatch objects
            string psmtsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MSMSids.psmtsv");
            TsvPsms = PsmTsvReader.ReadTsv(psmtsvPath, out var warnings);
            Psms = new List<PeptideSpectralMatch>();
            ProteinList = new List<Protein>();
            MyFileManager = new MyFileManager(true);

            foreach (PsmFromTsv readPsm in TsvPsms)
            {
                string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "TestData", "MbrAnalysisTest", readPsm.FileNameWithoutExtension + ".mzML");
                MsDataScan scan = MyFileManager.LoadFile(filePath, new CommonParameters()).GetOneBasedScan(readPsm.Ms2ScanNumber);
                Ms2ScanWithSpecificMass ms2Scan = new Ms2ScanWithSpecificMass(scan, readPsm.PrecursorMz, readPsm.PrecursorCharge,
                    filePath, new CommonParameters());
                Protein protein = new Protein(readPsm.BaseSeq, readPsm.ProteinAccession, readPsm.OrganismName,
                    isDecoy: readPsm.DecoyContamTarget == "D" ? true : false,
                    isContaminant: readPsm.DecoyContamTarget == "C" ? true : false);
                string[] startAndEndResidues = readPsm.StartAndEndResiduesInProtein.Split(" ");
                int startResidue = Int32.Parse(startAndEndResidues[0].Trim('['));
                int endResidue = Int32.Parse(startAndEndResidues[2].Trim(']'));

                PeptideWithSetModifications pwsm = new PeptideWithSetModifications(
                    readPsm.FullSequence, null, p: protein, digestionParams: new DigestionParams(),
                    oneBasedStartResidueInProtein: startResidue, oneBasedEndResidueInProtein: endResidue);
                PeptideSpectralMatch psm = new PeptideSpectralMatch(pwsm, 0, readPsm.Score, readPsm.Ms2ScanNumber, ms2Scan,
                    new CommonParameters(), readPsm.MatchedIons);
                psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                if (readPsm.Ms2ScanNumber == 206 && readPsm.BaseSeq.Equals("HADIVTTTTHK")) psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0.0046, 0); // Necessary for to be implemented "original pep" test
                Psms.Add(psm);
                ProteinList.Add(protein);
            }

            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput"));
        }

        private IEnumerable<IonQuantInfo> ParseCombinedIonsFile(string[] combinedIonLine)
        {
            IonQuantInfo ion = null;
            foreach (var kvp in SampleToFileDict)
            {
                int intensityIndex = CombinedIonHeader[kvp.Key + " Intensity"];
                int spectralCountIndex = CombinedIonHeader[kvp.Key + " Spectral Count"];
                if (Double.TryParse(combinedIonLine[intensityIndex], out double intensity) &&
                    Int32.TryParse(combinedIonLine[spectralCountIndex], out int spectralCount))
                {
                    if (intensity > 0 && spectralCount == 0)
                    {
                        try
                        {
                            ion = new IonQuantInfo(combinedIonLine, file: kvp.Value, intensity);
                        }
                        catch (ArgumentException e)
                        {
                            ion = null;
                        }
                        yield return (ion);
                    }
                }
            }

        }

        /// <summary>
        /// Reads in a list of IonQuantInfo objects and finds their corresponding psm file, then updates the retention time for each
        /// This is Hard Coded !!!
        /// However, it could be generalized with more string parsing
        /// </summary>
        /// <param name="ionQuantMbrPeptides"> List of IonQuantInfo objects generated from a combined_ion.tsv fragger file </param>
        private void UpdateIons(List<IonQuantInfo> ionQuantMbrPeptides)
        {
            // Psm info is stored in a file specific fashion

            string superFolder = @"D:\HelaSingleCellQCmzML\Fragger3";
            var ionsByFile = ionQuantMbrPeptides.OrderBy(i => i.BaseSequence).
                GroupBy(i => i.File);
            foreach (IGrouping<string, IonQuantInfo> ionGroup in ionsByFile)
            {
                // Fetch sample name from .mzML file name
                // In the future, sample names could be gleaned from the combined ion header and associated with each object
                // This would obviate the reverse dictionary lookup
                string sampleName = SampleToFileDict.Where(kvp => kvp.Value.Equals(ionGroup.Key)).Select(kvp => kvp.Key)
                    .First();
                string psmFilePath = Path.Join(superFolder, sampleName, "psm.tsv");
                IonQuantInfo[] ionArray = ionGroup.Skip(1).Select(i => i).ToArray();

                StreamReader reader = new StreamReader(psmFilePath);
                using (reader)
                {
                    // Remove header
                    reader.ReadLine();

                    while (!reader.EndOfStream)
                    {
                        string[] psmLine = reader.ReadLine().Split('\t');
                        string peptideSequence = psmLine[2] as string;
                        if (peptideSequence != null)
                        {
                            int ionIndex = Array.BinarySearch(ionArray, peptideSequence);
                            if (ionIndex >= 0)
                            {
                                ionArray[ionIndex].UpdateFileSpecificData(psmLine);
                            }
                        }
                    }
                }
            }
        }

    }

}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using iText.Kernel.Pdf;
using Microsoft.VisualBasic.CompilerServices;

namespace Test
{
    public class IonQuantInfo : IComparable
    {
        public string File { get; set; }
        public int Ms1ScanNumber { get; set; }
        public double PrecursorMz { get; set; }
        public double PrecursorMass { get; set; }
        public int PrecursorCharge { get; set; }
        public double RtMinutes { get; set; }
        public double Intensity { get; set; }
        public int? Ms2ScanNumber { get; set; }
        


        public string BaseSequence { get; set; }
        public string ModifiedSequence { get; set; }
        public string ProteinAccession { get; set; }
        public string Protein { get; set; }
        public int StartResidueInProtein { get; set; }
        public int EndResidueInProtein { get; set; }

        public IonQuantInfo(string file, int ms1ScanNumber, double mz, double precursorMass, int precursorCharge,
            double rTMinutes, double intensity, int? ms2ScanNumber, string baseSequence, string modifiedSequence, 
            string proteinAccession, string protein, int endResidueInProtein, int startResidueInProtein)
        {
            File = file;
            Ms1ScanNumber = ms1ScanNumber;
            PrecursorMz = mz;
            PrecursorMass = precursorMass;
            PrecursorCharge = precursorCharge;
            RtMinutes = rTMinutes;
            Intensity = intensity;
            Ms2ScanNumber = ms2ScanNumber;
            StartResidueInProtein = startResidueInProtein;
            BaseSequence = baseSequence;
            ModifiedSequence = modifiedSequence;
            ProteinAccession = proteinAccession;
            Protein = protein;
            EndResidueInProtein = endResidueInProtein;
        }

        public IonQuantInfo(string[] ionQuantLine, string file, double intensity)
        {
            if(!String.IsNullOrEmpty(ionQuantLine[0]) &&
                Int32.TryParse(ionQuantLine[4], out int startResidue) &&
                Int32.TryParse(ionQuantLine[5], out int endResidue) &&
                Double.TryParse(ionQuantLine[7], out double mz) &&
                Int32.TryParse(ionQuantLine[8], out int charge)
              )
            {
                BaseSequence = ionQuantLine[0];
                ModifiedSequence = ionQuantLine[1];
                ProteinAccession = ionQuantLine[11];
                Protein = ionQuantLine[12];
                StartResidueInProtein = startResidue;
                EndResidueInProtein = endResidue;
                PrecursorMz = mz;
                PrecursorCharge = charge;

                File = file;
                Intensity = intensity;
            }
            else
            {
                throw new ArgumentException();
            }

        }


        public int CompareTo(object? obj)
        {
            IonQuantInfo comparedIon = obj as IonQuantInfo;
            if (comparedIon != null && comparedIon.BaseSequence != null && this.BaseSequence != null)
            {
                return String.Compare(this.BaseSequence, comparedIon.BaseSequence);
            }

            String peptideSequence = obj as String;
            if (peptideSequence != null && this.BaseSequence != null)
            {
                return String.Compare(this.BaseSequence, peptideSequence);
            }

            throw new ArgumentException();
        }

        public void UpdateFileSpecificData(string[] psmLine)
        {
            if (Double.TryParse(psmLine[8], out double rtSeconds) &&
                Double.TryParse(psmLine[9], out double precursorMass) &&
                File.Equals(psmLine[1]))
            {
                RtMinutes = rtSeconds / 60.0;
                PrecursorMass = precursorMass;
            }
        }
    }
}

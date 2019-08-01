﻿using FlashLFQ;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class SilacConversions
    {
        private static readonly string LABEL_DELIMITER = " & ";

        public static string GetLabeledBaseSequence(string unlabeledBaseSequence, SilacLabel label)
        {
            if (label != null)
            {
                string labeledBaseSequence = unlabeledBaseSequence.Replace(label.OriginalAminoAcid, label.AminoAcidLabel);
                if (label.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                    {
                        labeledBaseSequence = labeledBaseSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel);
                    }
                }
                return labeledBaseSequence;
            }
            else
            {
                return unlabeledBaseSequence;
            }
        }

        public static PeptideSpectralMatch GetLabeledPsm(PeptideSpectralMatch psm, int notch, PeptideWithSetModifications pwsm, string labeledBaseSequence)
        {
            PeptideWithSetModifications labeledPwsm = new PeptideWithSetModifications(
                pwsm.Protein,
                pwsm.DigestionParams,
                pwsm.OneBasedStartResidueInProtein,
                pwsm.OneBasedEndResidueInProtein,
                pwsm.CleavageSpecificityForFdrCategory,
                pwsm.PeptideDescription,
                pwsm.MissedCleavages,
                pwsm.AllModsOneIsNterminus,
                pwsm.NumFixedMods,
                labeledBaseSequence);
            return psm.Clone(new List<(int Notch, PeptideWithSetModifications Peptide)> { (notch, labeledPwsm) });

        }

        public static PeptideSpectralMatch GetSilacPsm(PeptideSpectralMatch psm, SilacLabel silacLabel)
        {
            if (silacLabel == null)
            {
                return psm;
            }
            else
            {
                List<(int Notch, PeptideWithSetModifications Peptide)> updatedBestMatchingPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
                foreach ((int Notch, PeptideWithSetModifications Peptide) notchAndPwsm in psm.BestMatchingPeptides)
                {
                    PeptideWithSetModifications modifiedPwsm = CreateSilacPwsm(silacLabel, notchAndPwsm.Peptide);
                    updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
                }
                return psm.Clone(updatedBestMatchingPeptides);
            }
        }

        //modify the proteins to appear only light (we want a protein sequence to look like PROTEINK instead of PROTEINa)
        public static List<PeptideSpectralMatch> UpdateProteinSequencesToLight(List<PeptideSpectralMatch> originalPsms, List<SilacLabel> labels)
        {
            List<PeptideSpectralMatch> psmsToReturn = new List<PeptideSpectralMatch>();
            foreach(PeptideSpectralMatch psm in originalPsms)
            {
                List<(int Notch, PeptideWithSetModifications Peptide)> originalPeptides = psm.BestMatchingPeptides.ToList();
                List<(int Notch, PeptideWithSetModifications Peptide)> updatedPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
                foreach ((int Notch, PeptideWithSetModifications Peptide) notchPwsm in originalPeptides)
                {
                    PeptideWithSetModifications pwsm = notchPwsm.Peptide;
                    SilacLabel label = GetRelevantLabelFromBaseSequence(pwsm.BaseSequence, labels);
                    Protein updatedProtein = pwsm.Protein;
                    if (label != null)
                    {
                        string proteinLightSequence = updatedProtein.BaseSequence;
                        proteinLightSequence = proteinLightSequence.Replace(label.AminoAcidLabel, label.OriginalAminoAcid);
                        if (label.AdditionalLabels != null)
                        {
                            foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                            {
                                proteinLightSequence = proteinLightSequence.Replace(additionalLabel.AminoAcidLabel, additionalLabel.OriginalAminoAcid);
                            }
                        }
                        updatedProtein = new Protein(pwsm.Protein, proteinLightSequence);
                    }
                    PeptideWithSetModifications updatedPwsm = new PeptideWithSetModifications(
                        updatedProtein,
                        pwsm.DigestionParams,
                        pwsm.OneBasedStartResidueInProtein,
                        pwsm.OneBasedEndResidueInProtein,
                        pwsm.CleavageSpecificityForFdrCategory,
                        pwsm.PeptideDescription,
                        pwsm.MissedCleavages,
                        pwsm.AllModsOneIsNterminus,
                        pwsm.NumFixedMods,
                        pwsm.BaseSequence);
                    updatedPeptides.Add((notchPwsm.Notch, updatedPwsm));
                }

                psmsToReturn.Add(psm.Clone(updatedPeptides));
            }

            return psmsToReturn;
        }

        //Converts the heavy char label "a" into a human readable label "K+8.014"
        public static string GetSilacLightBaseSequence(string baseSequence, SilacLabel label)
        {
            if (label != null)
            {
                baseSequence = baseSequence.Replace(label.AminoAcidLabel.ToString(), HeavyStringForPeptides(label));
                if (label.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                    {
                        baseSequence = baseSequence.Replace(additionalLabel.AminoAcidLabel.ToString(), HeavyStringForPeptides(additionalLabel));
                    }
                }
            }
            return baseSequence;
        }

        //Converts the heavy char label "a" into a human readable label "K+8.014", but doesn't change mods containing the char like "Oxid'a'tion"
        public static string GetSilacLightFullSequence(string fullSequence, SilacLabel label, bool includeMassDifference = true)
        {
            //overwrite full sequence
            if (label != null)
            {
                List<SilacLabel> labels = new List<SilacLabel> { label };
                if (label.AdditionalLabels != null)
                {
                    labels.AddRange(label.AdditionalLabels);
                }

                foreach (SilacLabel additionalLabel in labels)
                {
                    string replacementSequence = includeMassDifference ? HeavyStringForPeptides(additionalLabel) : additionalLabel.OriginalAminoAcid.ToString();

                    bool inModification = false;
                    for (int i = 0; i < fullSequence.Length; i++)
                    {
                        if (inModification)
                        {
                            if (fullSequence[i] == ']')
                            {
                                inModification = false;
                            }
                        }
                        else
                        {
                            char currentChar = fullSequence[i];
                            if (currentChar == '[')
                            {
                                inModification = true;
                            }
                            else if (currentChar == additionalLabel.AminoAcidLabel)
                            {
                                fullSequence = fullSequence.Substring(0, i) + replacementSequence + fullSequence.Substring(i + 1, fullSequence.Length - i - 1);
                                i += replacementSequence.Length - 1; //-1 because we removed the label amino acid
                            }
                        }
                    }
                }
            }
            return fullSequence;
        }

        public static string GetAmbiguousLightSequence(string originalSequence, List<SilacLabel> labels, bool baseSequence)
        {
            string[] multipleSequences = originalSequence.Split('|').ToArray(); //if ambiguity
            string localSequence = "";
            foreach (string sequence in multipleSequences)
            {
                SilacLabel label = GetRelevantLabelFromBaseSequence(sequence, labels);
                localSequence += (baseSequence ? GetSilacLightBaseSequence(sequence, label) : GetSilacLightFullSequence(sequence, label)) + "|";
            }
            if (localSequence.Length != 0)
            {
                localSequence = localSequence.Substring(0, localSequence.Length - 1); //remove last "|"
            }
            return localSequence;
        }

        public static string HeavyStringForPeptides(SilacLabel label)
        {
            return label.OriginalAminoAcid + "(" + label.MassDifference + ")";
        }

        public static SilacLabel GetRelevantLabelFromBaseSequence(string baseSequence, List<SilacLabel> labels)
        {
            return labels.Where(x => baseSequence.Contains(x.AminoAcidLabel) ||
                (x.AdditionalLabels != null && x.AdditionalLabels.Any(y => baseSequence.Contains(y.AminoAcidLabel)))).FirstOrDefault();
        }

        public static SilacLabel GetRelevantLabelFromFullSequence(string fullSequence, List<SilacLabel> labels)
        {
            bool inModification = false;
            List<char> labelResidues = labels.Select(x => x.AminoAcidLabel).ToList();
            for (int i = 0; i < fullSequence.Length; i++)
            {
                if (inModification)
                {
                    if (fullSequence[i] == ']')
                    {
                        inModification = false;
                    }
                }
                else
                {
                    char currentChar = fullSequence[i];
                    if (currentChar == '[')
                    {
                        inModification = true;
                    }
                    else if (labelResidues.Contains(currentChar))
                    {
                        return labels.Where(x => currentChar == x.AminoAcidLabel).First();
                    }
                }
            }
            //if nothing
            return null;
        }

        public static PeptideWithSetModifications CreateSilacPwsm(SilacLabel silacLabel, PeptideWithSetModifications pwsm)
        {
            string baseSequence = pwsm.BaseSequence;
            Protein updatedProtein = null;

            baseSequence = baseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid); //create light sequence
            string proteinBaseSequence = pwsm.Protein.BaseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid);
            if (silacLabel.AdditionalLabels != null)
            {
                foreach (SilacLabel additionalLabel in silacLabel.AdditionalLabels)
                {
                    baseSequence = baseSequence.Replace(additionalLabel.AminoAcidLabel, additionalLabel.OriginalAminoAcid); //create light sequence
                    proteinBaseSequence = pwsm.Protein.BaseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid);
                }
            }

            return new PeptideWithSetModifications(
                updatedProtein ?? pwsm.Protein,
                pwsm.DigestionParams,
                pwsm.OneBasedStartResidueInProtein,
                pwsm.OneBasedEndResidueInProtein,
                pwsm.CleavageSpecificityForFdrCategory,
                pwsm.PeptideDescription,
                pwsm.MissedCleavages,
                pwsm.AllModsOneIsNterminus,
                pwsm.NumFixedMods,
                baseSequence);
        }

        public static SilacLabel AssignValidHeavyCharacter(SilacLabel originalLabel, char heavyLabel)
        {
            double massDifference = Convert.ToDouble(originalLabel.MassDifference.Substring(1));
            if (originalLabel.MassDifference[0] == '-')
            {
                massDifference *= -1;
            }
            //Add the silac residues to the dictionary
            Residue.AddNewResiduesToDictionary(new List<Residue> { new Residue(originalLabel.MassDifference, heavyLabel, heavyLabel.ToString(), Chemistry.ChemicalFormula.ParseFormula(originalLabel.LabelChemicalFormula), ModificationSites.All) });

            return new SilacLabel(originalLabel.OriginalAminoAcid, heavyLabel, originalLabel.LabelChemicalFormula, massDifference);
        }

        public static SpectraFileInfo GetHeavyFileInfo(SpectraFileInfo originalFile, SilacLabel label)
        {
            if (label == null)
            {
                return originalFile;
            }
            else
            {
                string heavyFileName = originalFile.FilenameWithoutExtension + "(" + label.OriginalAminoAcid + label.MassDifference;
                if (label.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionaLabel in label.AdditionalLabels)
                    {
                        heavyFileName += LABEL_DELIMITER + additionaLabel.OriginalAminoAcid + additionaLabel.MassDifference;
                    }
                }
                heavyFileName += ")." + originalFile.FullFilePathWithExtension.Split('.').Last(); //add extension

                return new SpectraFileInfo(heavyFileName, originalFile.Condition, originalFile.BiologicalReplicate, originalFile.TechnicalReplicate, originalFile.Fraction);
            }
        }

        public static HashSet<FlashLFQ.ProteinGroup> CleanPastProteinQuant(HashSet<FlashLFQ.ProteinGroup> originalProteinGroups)
        {
            HashSet<FlashLFQ.ProteinGroup> cleanedProteinGroups = new HashSet<FlashLFQ.ProteinGroup>();
            foreach(FlashLFQ.ProteinGroup pg in originalProteinGroups)
            {
                cleanedProteinGroups.Add(new FlashLFQ.ProteinGroup(pg.ProteinGroupName, pg.GeneName, pg.Organism));
            }
            return cleanedProteinGroups;
        }

        //If SILAC (Post-Quantification), compress the light/heavy protein group pairs into the same light protein group but different files
        //Create new files for each silac label and file so that "file 1" now becomes "file 1 (light)" and "file 1 (heavy)"
        //Change heavy residue into the light residue plus a string label ("PEPTIDEa" -> "PEPTIDEK(+8.014)")
        //This light to heavy conversion needs to happen for the flashLFQ peptides here, but can't for the psm peptides, which are constrained to the protein
        //i.e. pwsms currently don't have sequences; they have start/end residues and a protein sequence. We have to change the output sequences when they're created.
        public static void SilacConversionsPostQuantification(List<SilacLabel> allSilacLabels, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels,
            List<SpectraFileInfo> spectraFileInfo, List<ProteinGroup> proteinGroups, HashSet<DigestionParams> listOfDigestionParams, FlashLfqResults flashLfqResults,
            List<PeptideSpectralMatch> allPsms, Dictionary<string, int> modsToWriteSelection, bool quantifyUnlabeledPeptides)
        {
            //do protein quant if we had any results
            //if no results, we still may need to edit the psms
            if (flashLfqResults != null) //can be null if no unambiguous psms were found
            {
                //after this point, we now have quantification values for the peptides, but they all belong to the same "unlabeled" protein and are in the same file
                //We can remove "labeled" peptides from each file and put them in a new file as "unlabeled".

                //MAKE NEW RAW FILES
                //update number of spectra files to include a new file for each label*condition
                List<SpectraFileInfo> silacSpectraFileInfo = new List<SpectraFileInfo>(); //new files
                Dictionary<SpectraFileInfo, List<SpectraFileInfo>> originalToLabeledFileInfoDictionary = new Dictionary<SpectraFileInfo, List<SpectraFileInfo>>();
                List<SpectraFileInfo> unlabeledFiles = new List<SpectraFileInfo>(); //null can't be a key in the dictionary, so this list is used for the unlabeled files.

                flashLfqResults.SpectraFiles.Clear(); //clear existing files so we can replace them with labeled ones

                //foreach existing file
                //if multiplex
                if (turnoverLabels == null)
                {
                    //populate dictionary
                    if (quantifyUnlabeledPeptides)
                    {
                        spectraFileInfo.ForEach(x => originalToLabeledFileInfoDictionary.Add(x, new List<SpectraFileInfo> { x }));
                        //get the unlabeled
                        unlabeledFiles = spectraFileInfo;
                        flashLfqResults.SpectraFiles.AddRange(spectraFileInfo);
                    }
                    else
                    {
                        spectraFileInfo.ForEach(x => originalToLabeledFileInfoDictionary.Add(x, new List<SpectraFileInfo>()));
                    }

                    //get the labeled
                    foreach (SilacLabel label in allSilacLabels)
                    {
                        List<SpectraFileInfo> labeledFiles = new List<SpectraFileInfo>();
                        foreach (SpectraFileInfo originalFile in spectraFileInfo)
                        {
                            //foreach label, add a new file with the label                      
                            //silacSpectraFileInfo.Add(silacFile);
                            SpectraFileInfo labeledInfo = SilacConversions.GetHeavyFileInfo(originalFile, label);
                            labeledFiles.Add(labeledInfo);
                            originalToLabeledFileInfoDictionary[originalFile].Add(labeledInfo);
                        }
                        flashLfqResults.SpectraFiles.AddRange(labeledFiles);
                    }
                }

                //we have the files, now let's reassign the psms.
                //there are a few ways to do this, but we're going to generate the "base" peptide and assign to that

                //Dictionary of protein accession to peptides
                Dictionary<string, List<(SilacLabel label, FlashLFQ.Peptide peptide)>> unlabeledToPeptidesDictionary = new Dictionary<string, List<(SilacLabel label, FlashLFQ.Peptide peptide)>>();
                Dictionary<string, FlashLFQ.Peptide> updatedPeptideModifiedSequences = new Dictionary<string, FlashLFQ.Peptide>();
                IEnumerable<FlashLFQ.Peptide> lfqPwsms = flashLfqResults.PeptideModifiedSequences.Values;
                if (turnoverLabels == null || turnoverLabels.Value.startLabel == null || turnoverLabels.Value.endLabel == null) //if multiplex, or if at least one of the turnover labels is unlabeled
                {
                    foreach (FlashLFQ.Peptide peptide in lfqPwsms)
                    {
                        //convert to the unlabeled sequence
                        string labeledSequence = peptide.Sequence;
                        SilacLabel label = GetRelevantLabelFromFullSequence(labeledSequence, allSilacLabels);
                        string unlabeledSequence = GetSilacLightFullSequence(labeledSequence, label, false);
                        if (unlabeledToPeptidesDictionary.ContainsKey(unlabeledSequence))
                        {
                            unlabeledToPeptidesDictionary[unlabeledSequence].Add((label, peptide));
                        }
                        else
                        {
                            unlabeledToPeptidesDictionary.Add(unlabeledSequence, new List<(SilacLabel, FlashLFQ.Peptide)> { (label, peptide) });
                        }
                    }
                }
                else //if both start and end labels exist for turnover, then it can be trickier to find the start label
                {
                    SilacLabel startLabel = turnoverLabels.Value.startLabel;
                    SilacLabel endLabel = turnoverLabels.Value.endLabel;
                    foreach (FlashLFQ.Peptide peptide in lfqPwsms)
                    {
                        string originalSequence = peptide.Sequence;
                        string partiallyCleanedSequence = originalSequence;
                        SilacLabel labelToReport = null;
                        //convert to the unlabeled sequence
                        if (endLabel != null)
                        {
                            partiallyCleanedSequence = SilacConversions.GetSilacLightFullSequence(originalSequence, endLabel, false);
                            if (!partiallyCleanedSequence.Equals(originalSequence)) //if any residue is new, then the whole peptide must be newly synthesized
                            {
                                labelToReport = endLabel;
                            }
                        }
                        string fullyCleanedSequence = partiallyCleanedSequence;
                        if (startLabel != null) //there might also be some old residues in the new peptide, so this is an "if" and not an "else"
                        {
                            fullyCleanedSequence = SilacConversions.GetSilacLightFullSequence(partiallyCleanedSequence, startLabel, false);
                            if (!fullyCleanedSequence.Equals(partiallyCleanedSequence) && labelToReport == null)
                            {
                                labelToReport = startLabel;
                            }
                        }

                        if (unlabeledToPeptidesDictionary.ContainsKey(fullyCleanedSequence))
                        {
                            unlabeledToPeptidesDictionary[fullyCleanedSequence].Add((labelToReport, peptide));
                        }
                        else
                        {
                            unlabeledToPeptidesDictionary.Add(fullyCleanedSequence, new List<(SilacLabel, FlashLFQ.Peptide)> { (labelToReport, peptide) });
                        }
                    }
                }

                //we now have a dictionary of unlabeledBaseSequence to the labeled peptides
                //Better SILAC results can be obtained by using the summed intensities from ms1 scans where all peaks were found, rather than the apex
                //foreach peptide, unlabeled peptide, get the isotopic envelope intensities for each labeled peptide in each file
                //save the intensities from ms1s that are shared. If no ms1s contains all the peaks, then just use the apex intensity (default)
                foreach (KeyValuePair<SpectraFileInfo, List<ChromatographicPeak>> kvp in flashLfqResults.Peaks)
                {
                    //make a dictionary for easy peak lookup from peptide sequences
                    Dictionary<string, ChromatographicPeak> sequenceToPeakDictionary = new Dictionary<string, ChromatographicPeak>();
                    foreach (ChromatographicPeak peak in kvp.Value)
                    {
                        sequenceToPeakDictionary.Add(peak.Identifications.First().ModifiedSequence, peak);
                    }

                    foreach (List<(SilacLabel label, FlashLFQ.Peptide peptide)> peptideGroup in unlabeledToPeptidesDictionary.Select(x => x.Value))
                    {
                        List<string> sequences = peptideGroup.Select(x => x.peptide.Sequence).ToList();

                        //get peaks of interest for this peptide group
                        List<ChromatographicPeak> peaks = kvp.Value;
                        List<ChromatographicPeak> peaksOfInterest = new List<ChromatographicPeak>();
                        foreach (string sequence in sequences)
                        {
                            if (sequenceToPeakDictionary.TryGetValue(sequence, out ChromatographicPeak peakOfInterest))
                            {
                                peaksOfInterest.Add(peakOfInterest);
                            }
                        }

                        //If there are fewer than 2 peaks, we can't do any comparisons
                        if (peaksOfInterest.Count > 1)
                        {
                            //get isotopic envelopes that are shared by all
                            List<int> scanIndex = peaksOfInterest.First().IsotopicEnvelopes.Select(x => x.IndexedPeak).Select(x => x.ZeroBasedMs1ScanIndex).ToList();
                            for (int i = 1; i < peaksOfInterest.Count; i++)
                            {
                                List<int> currentScanIndexes = peaksOfInterest[i].IsotopicEnvelopes.Select(x => x.IndexedPeak).Select(x => x.ZeroBasedMs1ScanIndex).ToList();
                                scanIndex = scanIndex.Intersect(currentScanIndexes).ToList();
                                if (scanIndex.Count == 0) //if there's no overlap, then we're done!
                                {
                                    break;
                                }
                            }
                            //if we aren't sticking with the default values
                            if (scanIndex.Count != 0)
                            {
                                //update peptides
                                foreach (FlashLFQ.Peptide peptide in peptideGroup.Select(x => x.peptide))
                                {
                                    ChromatographicPeak peakForThisPeptide = peaksOfInterest.Where(x => peptide.Sequence.Equals(x.Identifications.First().ModifiedSequence)).First();
                                    double summedIntensity = peakForThisPeptide.IsotopicEnvelopes.Where(x => scanIndex.Contains(x.IndexedPeak.ZeroBasedMs1ScanIndex)).Select(x => x.Intensity).Sum();
                                    peptide.SetIntensity(kvp.Key, summedIntensity);
                                }
                            }
                        }
                    }
                }


                //SPLIT THE FILES
                List<FlashLFQ.Peptide> updatedPeptides = new List<FlashLFQ.Peptide>();

                //split the heavy/light peptides into separate raw files, remove the heavy peptide
                if (turnoverLabels != null)
                {
                    Dictionary<SpectraFileInfo, double[]> fileToRecycleDictionary = new Dictionary<SpectraFileInfo, double[]>();
                    //TODO just want the original files, no reason to look at the silac ones
                    foreach (SpectraFileInfo info in spectraFileInfo)
                    {
                        fileToRecycleDictionary[info] = new double[3];
                    }

                    //go through each peptide that has more than one label
                    foreach (KeyValuePair<string, List<(SilacLabel label, FlashLFQ.Peptide peptide)>> kvp in unlabeledToPeptidesDictionary.Where(x => x.Value.Count > 2))
                    {
                        //the order should always be from start to end
                        //TODO if that's not true, use the silac label
                        //TODO if that is true, remove the silaclabel from the dictionary
                        //TODO check that values are present for light and for double the label mass to prevent incorrectly id'd missed cleavages from swaying data?
                        Dictionary<SpectraFileInfo, List<double>> tempFileToRecycleDictionary = new Dictionary<SpectraFileInfo, List<double>>();
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            tempFileToRecycleDictionary[info] = new List<double>();
                        }

                        List<FlashLFQ.Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                        foreach (FlashLFQ.Peptide peptide in peptides) //assumes sorted old->new
                        {
                            foreach (SpectraFileInfo info in spectraFileInfo)
                            {
                                tempFileToRecycleDictionary[info].Add(peptide.GetIntensity(info));
                            }
                        }

                        //aggregate into LL, HL, and HH if they contain 2+ missed cleavages
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            List<double> values = tempFileToRecycleDictionary[info];
                            if (values.Count > 6)
                            {
                                //we can't really use data with this many missed cleavages (data isn't super accurate), so remove them
                                values = new List<double> { 0, 0, 0 };
                            }
                            else if (values.Count > 3)
                            {
                                //TODO: make sure this works with a ton of missed cleavages (5)
                                if (values.Count == 6)
                                {
                                    //LLLLL, LLLLH, LLLHH, LLHHH, LHHHH, HHHHHH have n combinations of 1, 5, 10, 10, 5, 1
                                    double LLLL = values[0] + values[1] / 5;
                                    double LLLH = values[1] * 4 / 5 + values[2] * 4 / 10;
                                    double LLHH = values[2] * 6 / 10 + values[3] * 6 / 10;
                                    double LHHH = values[3] * 4 / 10 + values[4] * 4 / 5;
                                    double HHHH = values[4] / 5 + values[5];
                                    values = new List<double> { LLLL, LLLH, LLHH, LHHH, HHHH };
                                }
                                if (values.Count == 5)
                                {
                                    //LLLL, LLLH, LLHH, LHHH, HHHH have n combinations of 1, 4, 6, 4, 1
                                    double LLL = values[0] + values[1] / 4;
                                    double LLH = values[1] * 3 / 4 + values[2] / 2; // divide by two is multiply by three divide by six
                                    double LHH = values[2] / 2 + values[3] * 3 / 4;
                                    double HHH = values[3] / 4 + values[4];
                                    values = new List<double> { LLL, LLH, LHH, HHH };
                                }
                                if (values.Count == 4)
                                {
                                    //LLL, LLH, LHH, HHH have probabilities of 1, 3, 3, 1
                                    double LL = values[0] + values[1] / 3; // divide by three for the probability difference (only one way to make LL from LLH)
                                    double LH = values[1] * 2 / 3 + values[2] * 2 / 3; //divide by three, multiply by two (two ways to make LH from LLH or from LHH)
                                    double HH = values[2] / 3 + values[3]; //same logic as LL
                                    values = new List<double> { LL, LH, HH };
                                }

                                //TODO: validate this
                                //overwrite the current peptide quantification values with the total values shown here IF there was more than one missed cleavage
                                for (int i = 0; i < 3; i++)
                                {
                                    peptides[i].SetIntensity(info, values[i]);
                                    //TODO: does the detection type stay here, or do we need to set it?
                                }
                            }
                            //else if three values, then we're fine
                            //add the values from this peptide to the file-specific-values of all peptides
                            double[] totalValues = fileToRecycleDictionary[info];
                            for (int i = 0; i < values.Count; i++)
                            {
                                totalValues[i] += values[i];
                            }
                        }
                    }

                    //with the summed LL, LH, and HH values for each file, we can calculate the probability that a heavy amino acid is incorporated (Ph)
                    //the probability a light amino acid is incorporated (Pl) is simply Pl = 1 - Ph
                    //We can simplify the math if we assume that all protein started as completely light, @t=0: LL=100%, LH = 0%, and HH = 0%
                    //N = 1 - Qll + Qlh * Qlh / (4 * Qhh), (where N is the fraction of all protein that has turned over and Qll is the Quantity(fraction) of LL at the current timepoint
                    //Qhh = qhh -N*qhh+Ph*Ph*N, where Q is quantity(fraction) at the timepoint and q is the quantity(fraction) at the start point (time zero).
                    //Ph = sqrt((Qhh-qhh+N*qhh)/N), and qhh = 0, so we can simplify as
                    //Ph = sqrt(Qhh/N)
                    Dictionary<SpectraFileInfo, double> fileToHeavyProbabilityDictionary = new Dictionary<SpectraFileInfo, double>();
                    foreach (SpectraFileInfo info in spectraFileInfo)
                    {
                        //convert the absolute values of LL, LH, and HH to relative values
                        double[] values = fileToRecycleDictionary[info];
                        double sum = values.Sum();
                        for (int i = 0; i < values.Length; i++)
                        {
                            values[i] = values[i] / sum;
                        }
                        double Ph = Math.Sqrt(values[2] / (1 - values[0] + Math.Pow(values[1], 2) / (4 * values[2])));
                        fileToHeavyProbabilityDictionary[info] = Ph;
                    }


                    //foreach group, the labeled peptides should be split into their labeled files
                    //we're deleting the heavy results after we pull those results into a different file
                    //additionally, the Ph will be used to correct both values (L is originally too high, H too low)
                    foreach (SpectraFileInfo info in spectraFileInfo)
                    {
                        string fullPathWithExtension = info.FullFilePathWithExtension;
                        string[] pathArray = fullPathWithExtension.Split('.');
                        string extension = pathArray.Last();
                        string filePathWithoutExtension = fullPathWithExtension.Substring(0, fullPathWithExtension.Length - extension.Length - 1); //-1 removes the '.'
                        SpectraFileInfo lightInfo = new SpectraFileInfo(filePathWithoutExtension + "_Old." + extension, info.Condition, info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                        SpectraFileInfo heavyInfo = new SpectraFileInfo(filePathWithoutExtension + "_New." + extension, info.Condition, info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                        originalToLabeledFileInfoDictionary[info] = new List<SpectraFileInfo> { lightInfo, heavyInfo };
                        flashLfqResults.SpectraFiles.Add(lightInfo);
                        flashLfqResults.SpectraFiles.Add(heavyInfo);
                    }
                    foreach (KeyValuePair<string, List<(SilacLabel label, FlashLFQ.Peptide peptide)>> kvp in unlabeledToPeptidesDictionary)
                    {
                        string unlabeledSequence = kvp.Key; //this will be the key for the new quant entry
                        List<FlashLFQ.Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                        //if no missed cleavages
                        FlashLFQ.Peptide lightPeptide;
                        FlashLFQ.Peptide mixedPeptide = null; //null assignment needed to build
                        FlashLFQ.Peptide heavyPeptide;
                        bool mixedPeptideExists = peptides.Count != 2;
                        if (!mixedPeptideExists)
                        {
                            lightPeptide = peptides[0];
                            heavyPeptide = peptides[1];
                        }
                        else
                        {
                            lightPeptide = peptides[0];
                            mixedPeptide = peptides[1];
                            heavyPeptide = peptides[2];
                        }

                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            double Ph = fileToHeavyProbabilityDictionary[info];
                            double light = lightPeptide.GetIntensity(info);
                            double heavy = heavyPeptide.GetIntensity(info);
                            double mixed = mixedPeptideExists ? mixedPeptide.GetIntensity(info) : 0;

                            //all the heavy is new, but some of the light is also new protein
                            //Ph helps out here. The correction factor is Pl*Qh/Ph, or (1-Ph)*Qh/Ph.
                            //it's possible that we obtain a negative value for the start, which doesn't make sense.
                            //This occurs when Qh>Ph. In these cases, set the light to zero and the heavy to light+heavy.
                            double correction = (1 - Ph) * heavy / Ph;
                            var updatedInfo = originalToLabeledFileInfoDictionary[info];
                            SpectraFileInfo startInfo = updatedInfo[0];
                            SpectraFileInfo endInfo = updatedInfo[1];

                            FlashLFQ.Peptide updatedPeptide = new FlashLFQ.Peptide(lightPeptide.Sequence, lightPeptide.UseForProteinQuant);
                            updatedPeptide.ProteinGroups = CleanPastProteinQuant(lightPeptide.ProteinGroups); //needed to keep protein info.

                            updatedPeptide.SetIntensity(startInfo, light - correction); //assign the corrected light intensity
                            updatedPeptide.SetDetectionType(startInfo, lightPeptide.GetDetectionType(info));
                            updatedPeptide.SetIntensity(endInfo, heavy + mixed + correction); //assign the corrected heavy intensity to the heavy file
                            updatedPeptide.SetDetectionType(endInfo, heavyPeptide.GetDetectionType(info)); //could include the mixed here if it really matters
                            updatedPeptides.Add(updatedPeptide);
                        }
                    }
                }
                else //multiplex
                {
                    //foreach original file
                    foreach (SpectraFileInfo info in spectraFileInfo)
                    {
                        foreach (var kvp in unlabeledToPeptidesDictionary)
                        {
                            string unlabeledSequence = kvp.Key;
                            List<FlashLFQ.Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                            FlashLFQ.Peptide representativePeptide = peptides[0];
                            FlashLFQ.Peptide updatedPeptide = new FlashLFQ.Peptide(unlabeledSequence, representativePeptide.UseForProteinQuant);
                            updatedPeptide.ProteinGroups = CleanPastProteinQuant(representativePeptide.ProteinGroups); //needed to keep protein info.
                            List<SpectraFileInfo> filesForThisFile = originalToLabeledFileInfoDictionary[info];
                            for (int i = 0; i < peptides.Count; i++) //the files and the peptides can use the same index, because there should be a distinct file for each label/peptide
                            {
                                SpectraFileInfo currentInfo = filesForThisFile[i];
                                FlashLFQ.Peptide currentPeptide = peptides[i];
                                updatedPeptide.SetIntensity(currentInfo, currentPeptide.GetIntensity(info));
                                updatedPeptide.SetDetectionType(currentInfo, currentPeptide.GetDetectionType(info));
                            }
                            updatedPeptides.Add(updatedPeptide);
                        }
                    }
                }
                //Update peptides
                var peptideResults = flashLfqResults.PeptideModifiedSequences;
                peptideResults.Clear();
                foreach (FlashLFQ.Peptide peptide in updatedPeptides)
                {
                    peptideResults.Add(peptide.Sequence, peptide);
                }

                //Do protein quant
                flashLfqResults.ProteinGroups.Clear();
                flashLfqResults.CalculateProteinResultsTop3();

                //update proteingroups to have all files for quantification
                if (proteinGroups != null)
                {
                    List<SpectraFileInfo> allInfo = originalToLabeledFileInfoDictionary.SelectMany(x => x.Value).ToList();
                    foreach (ProteinGroup proteinGroup in proteinGroups)
                    {
                        proteinGroup.FilesForQuantification = allInfo;
                        proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                        foreach (var spectraFile in proteinGroup.FilesForQuantification)
                        {
                            if (flashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
                            {
                                proteinGroup.IntensitiesByFile.Add(spectraFile, flashLfqProteinGroup.GetIntensity(spectraFile));
                            }
                            else
                            {
                                proteinGroup.IntensitiesByFile.Add(spectraFile, 0);
                            }
                        }
                    }
                }

                //Convert all lfqpeaks from heavy (a) to light (K+8.014) for output

                if (flashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                {
                    var lfqPeaks = flashLfqResults.Peaks;
                    List<SpectraFileInfo> peakKeys = lfqPeaks.Keys.ToList();

                    foreach (SpectraFileInfo key in peakKeys)
                    {
                        List<ChromatographicPeak> peaks = lfqPeaks[key];
                        for (int i = 0; i < peaks.Count; i++)
                        {
                            var peak = peaks[i];
                            List<Identification> identifications = new List<Identification>();
                            //check if we're removing light peaks and if it's a light peak
                            if (peak.Identifications.Any(x => SilacConversions.GetRelevantLabelFromBaseSequence(x.BaseSequence, allSilacLabels) != null)) //if no ids have any labels, remove them
                            {
                                List<Identification> updatedIds = new List<Identification>();
                                foreach (var id in peak.Identifications)
                                {
                                    SilacLabel label = SilacConversions.GetRelevantLabelFromBaseSequence(id.BaseSequence, allSilacLabels);

                                    Identification updatedId = new Identification(
                                        id.FileInfo,
                                        SilacConversions.GetSilacLightBaseSequence(id.BaseSequence, label),
                                        SilacConversions.GetSilacLightFullSequence(id.ModifiedSequence, label),
                                        id.MonoisotopicMass,
                                        id.Ms2RetentionTimeInMinutes,
                                        id.PrecursorChargeState,
                                        id.ProteinGroups.ToList(),
                                        id.OptionalChemicalFormula,
                                        id.UseForProteinQuant
                                        );
                                    updatedIds.Add(updatedId);
                                }
                                peak.Identifications.Clear();
                                peak.Identifications.AddRange(updatedIds);
                            }
                        }
                    }
                }
            }

            //convert all peaks into human readable format (PEPTIDE(K+8) and not PEPTIDEa)
            //convert all psms into human readable format
            //We can do this for all of the FlashLFQ peptides/peaks, because they use string sequences.
            //We are unable to do this for Parameters.AllPsms, because they store proteins and start/end residues instead
            //for Psms, we need to convert during the writing.
            for (int i = 0; i < allPsms.Count; i++)
            {
                allPsms[i].ResolveHeavySilacLabel(allSilacLabels, modsToWriteSelection);
            }
        }
    }
}
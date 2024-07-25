﻿using Chemistry;
using EngineLayer.CrosslinkSearch;
using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using Microsoft.ML.Data;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public static class PEP_Analysis_Cross_Validation
    {
        private static readonly double AbsoluteProbabilityThatDistinguishesPeptides = 0.05;
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

        // These are variables for the GBDT model
        private static readonly int _numberOfTrees = 1500; // Microsoft default is 100
        private static readonly double _learningRate = 0.2; // Microsoft default is 0.2
        public static double FdrThresholdForPositiveExamples { get; set; } = 0.005;

        public static string ComputePEPValuesForAllPSMsGeneric(List<SpectralMatch> psms, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, string outputFolder, FdrAnalysisEngine fdrEngine = null)
        {
            string[] trainingVariables = PsmData.trainingInfos[searchType];

            var qValueArray = psms.OrderByDescending(p => p).Select(p => p.FdrInfo.QValue).ToArray();
            FdrThresholdForPositiveExamples = qValueArray[(int)((double)qValueArray.Length * 0.25)];

            #region Construct Peptide Groups and Score Dictionaries

            List<PeptideMatchGroup> peptideGroups = new();
            foreach (var pepGroup in psms.Where(psm => psm.FullSequence != null).GroupBy(psm => psm.FullSequence))
            { 
                List<SpectralMatch> spectralMatches = pepGroup.ToList();
                peptideGroups.Add(new PeptideMatchGroup(pepGroup.Key, spectralMatches));
            }
            foreach(var psm in psms.Where(psm => psm.FullSequence == null))
            {
                peptideGroups.Add(new PeptideMatchGroup(psm.BaseSequence, new List<SpectralMatch> { psm }));
            }

            //ensure that the order is always stable.
            peptideGroups = peptideGroups
                .OrderByDescending(p => p.TargetCount)
                .ThenByDescending(p => p.DecoyCount)
                .OrderByDescending(p => p.BestMatch).ToList();

            int numberOfFilesRepresentedInPeptideList = peptideGroups.GroupBy(b => b.BestMatch.FullFilePath).ToList().Count();
            bool allFilesContainPeptides = numberOfFilesRepresentedInPeptideList == fileSpecificParameters.Count; //rare condition where each file has psms but some files don't have peptides. probably only happens in unit tests.

            bool trainOnPeptides = peptideGroups.Count() > 100;// && allFilesContainPeptides;
            if (!trainOnPeptides)
            {
                // If we have to train on the PSM level, then we need to reconstruct the list of "peptideGroups"
                // such that each psm is represented by its own peptide group
                peptideGroups = new();
                foreach (var psm in psms)
                {
                    peptideGroups.Add(new PeptideMatchGroup(psm.FullSequence ?? psm.BaseSequence, new List<SpectralMatch> { psm }));
                }
                peptideGroups = peptideGroups.OrderByDescending(p => p.BestMatch).ToList();
            }

            int chargeStateMode = chargeStateMode = GetChargeStateMode(peptideGroups.SelectMany(p => p));
            Dictionary<string, float> fileSpecificMedianFragmentMassErrors = GetFileSpecificMedianFragmentMassError(
                    peptideGroups.SelectMany(p => p));


            //These two dictionaries contain the average and standard deviations of hydrophobicitys measured in 1 minute increments accross each raw
            //file separately. An individully measured hydrobophicty calculated for a specific PSM sequence is compared to these values by computing
            //the z-score. That z-score is used as a feature for machine learning.
            //Separate dictionaries are created for peptides with modifications because SSRcalc doesn't really do a good job predicting hyrophobicity

            //The first string in the dictionary is the filename
            //The value of the dictionary is another dictionary that profiles the hydrophobicity behavior.
            //Each key is a retention time rounded to the nearest minute.
            //The value Tuple is the average and standard deviation, respectively, of the predicted hydrophobicities of the observed peptides eluting at that rounded retention time.
            if (trainingVariables.Contains("HydrophobicityZScore"))
            {
                fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = ComputeHydrophobicityValues(peptideGroups.SelectMany(p => p), fileSpecificParameters, false);
                fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = ComputeHydrophobicityValues(peptideGroups.SelectMany(p => p), fileSpecificParameters, true);
                fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = ComputeMobilityValues(peptideGroups.SelectMany(p => p), fileSpecificParameters);
            }

            // Partition the peptides into groups for cross-validation
            // the number of groups used for cross-validation is hard-coded at four. Do not change this number without changes other areas of effected code.
            int numGroups = 3;

            // Partition the peptide groups for cross-validation
            List<int>[] peptideGroupIndices = GetPeptideGroupIndices(peptideGroups, numGroups);

            #endregion Construct Peptide Groups and Score Dictionaries

            #region Construct Training Data

            //the psms will be randomly divided. but then we want to make another array that just contains the subset of peptides that are in those psms. that way we don't compute pep using any peptides that were used in training.
            IEnumerable<PsmData>[] PSMDataGroups = new IEnumerable<PsmData>[numGroups];

            int decoyCountNow = peptideGroups.SelectMany(p => p).Count(p => p.IsDecoy);

            int peptideGroupCount = peptideGroups.Count();
            int psmCounts = peptideGroups.SelectMany(p => p).Count();

            bool allSetsContainPositiveAndNegativeTrainingExamples = true;
            for (int i = 0; i < numGroups; i++)
            {
                PSMDataGroups[i] = CreatePsmData(searchType, fileSpecificParameters, peptideGroups, peptideGroupIndices[i], fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode);

                int targetCount = PSMDataGroups[i].Count(p => p.Label == true);
                int decoyCount = PSMDataGroups[i].Count(p => p.Label == false);

                if (!PSMDataGroups[i].Any(p => p.Label == true) || !PSMDataGroups[i].Any(p => p.Label == false))
                {
                    allSetsContainPositiveAndNegativeTrainingExamples = false;
                    break;
                }
            }

            MLContext mlContext = new MLContext();
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", 
                numberOfTrees: _numberOfTrees, learningRate: _learningRate);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();
            int sumOfAllAmbiguousPeptidesResolved = 0;

            #endregion

            #region Cross Validation/Scoring

            if (!allSetsContainPositiveAndNegativeTrainingExamples)
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";

            }
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;


            for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
            {
                List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                allGroupIndexes.RemoveAt(groupIndexNumber);

                //concat doesn't work in a loop, therefore I had to hard code the concat to group 2 out of 3 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                IDataView dataView = mlContext.Data.LoadFromEnumerable(PSMDataGroups[allGroupIndexes[0]]);
                if (numGroups > 2)
                {
                    dataView = mlContext.Data.LoadFromEnumerable(
                        PSMDataGroups[allGroupIndexes[0]]
                        .Concat(PSMDataGroups[allGroupIndexes[1]]));
                }
                trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(PSMDataGroups[groupIndexNumber]));
                CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                //Parallel operation of the following code requires the method to be stored and then read, once for each thread
                //if not output directory is specified, the model cannot be stored, and we must force single-threaded operation
                if (outputFolder != null)
                {
                    mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                }

                //model is trained on peptides but here we can use that to compute PEP for all PSMs
                int ambiguousPeptidesResolved = Compute_PSM_PEP(peptideGroups, peptideGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], searchType, fileSpecificParameters, fileSpecificMedianFragmentMassErrors, chargeStateMode, outputFolder);

                allMetrics.Add(metrics);
                sumOfAllAmbiguousPeptidesResolved += ambiguousPeptidesResolved;
            }

            #endregion

            

            return AggregateMetricsForOutput(allMetrics, sumOfAllAmbiguousPeptidesResolved);
        }

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics, int sumOfAllAmbiguousPeptidesResolved)
        {
            List<double> accuracy = allMetrics.Select(m => m.Accuracy).ToList();
            List<double> areaUnderRocCurve = allMetrics.Select(m => m.AreaUnderRocCurve).ToList();
            List<double> areaUnderPrecisionRecallCurve = allMetrics.Select(m => m.AreaUnderPrecisionRecallCurve).ToList();
            List<double> F1Score = allMetrics.Select(m => m.F1Score).ToList();
            List<double> logLoss = allMetrics.Select(m => m.LogLoss).ToList();
            List<double> logLossReduction = allMetrics.Select(m => m.LogLossReduction).ToList();
            List<double> positivePrecision = allMetrics.Select(m => m.PositivePrecision).ToList();
            List<double> positiveRecall = allMetrics.Select(m => m.PositiveRecall).ToList();
            List<double> negativePrecision = allMetrics.Select(m => m.NegativePrecision).ToList();
            List<double> negativeRecall = allMetrics.Select(m => m.NegativeRecall).ToList();

            // log-loss can stochastically take on a value of infinity.
            // correspondingly, log-loss reduction can be negative infinity.
            // when this happens for one or more of the metrics, it can lead to uninformative numbers.
            // so, unless they are all infinite, we remove them from the average. If they are all infinite, we report that.

            logLoss.RemoveAll(x => x == Double.PositiveInfinity);
            logLossReduction.RemoveAll(x => x == Double.NegativeInfinity);

            double logLossAverage = Double.PositiveInfinity;
            double logLossReductionAverage = Double.NegativeInfinity;

            if ((logLoss != null) && (logLoss.Any()))
            {
                logLossAverage = logLoss.Average();
            }

            if ((logLossReduction != null) && (logLossReduction.Any()))
            {
                logLossReductionAverage = logLossReduction.Average();
            }

            StringBuilder s = new StringBuilder();
            s.AppendLine();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + accuracy.Average().ToString());
            s.AppendLine("*       Area Under Curve:  " + areaUnderRocCurve.Average().ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + areaUnderPrecisionRecallCurve.Average().ToString());
            s.AppendLine("*       F1Score:  " + F1Score.Average().ToString());
            s.AppendLine("*       LogLoss:  " + logLossAverage.ToString());
            s.AppendLine("*       LogLossReduction:  " + logLossReductionAverage.ToString());
            s.AppendLine("*       PositivePrecision:  " + positivePrecision.Average().ToString());
            s.AppendLine("*       PositiveRecall:  " + positiveRecall.Average().ToString());
            s.AppendLine("*       NegativePrecision:  " + negativePrecision.Average().ToString());
            s.AppendLine("*       NegativeRecall:  " + negativeRecall.Average().ToString());
            s.AppendLine("*       Count of Ambiguous Peptides Removed:  " + sumOfAllAmbiguousPeptidesResolved.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        public static int Compute_PSM_PEP(List<PeptideMatchGroup> peptideGroups, List<int> peptideGroupIndices, MLContext mLContext, TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, string outputFolder)
        {
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;
            object lockObject = new object();
            int ambiguousPeptidesResolved = 0;

            //the trained model is not threadsafe. Therefore, to use the same model for each thread saved the model to disk. Then each thread reads its own copy of the model back from disk.
            //If there is no output folder specified, then this can't happen. We set maxthreads eqaul to one and use the model that gets passed into the method.
            if (String.IsNullOrEmpty(outputFolder))
            {
                maxThreads = 1;
            }

            Parallel.ForEach(Partitioner.Create(0, peptideGroupIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    ITransformer threadSpecificTrainedModel;
                    if (maxThreads == 1)
                    {
                        threadSpecificTrainedModel = trainedModel;
                    }
                    else
                    {
                        threadSpecificTrainedModel = mLContext.Model.Load(Path.Combine(outputFolder, "model.zip"), out DataViewSchema savedModelSchema);
                    }

                    // one prediction engine per thread, because the prediction engine is not thread-safe
                    var threadPredictionEngine = mLContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(threadSpecificTrainedModel);

                    int ambigousPeptidesRemovedinThread = 0;

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        foreach (SpectralMatch psm in peptideGroups[peptideGroupIndices[i]])
                        {
                            if (psm == null) continue;

                            List<int> indiciesOfPeptidesToRemove = new List<int>();
                            List<double> pepValuePredictions = new List<double>();

                            //Here we compute the pepvalue predection for each ambiguous peptide in a PSM. Ambiguous peptides with lower pepvalue predictions are removed from the PSM.

                            List<int> allBmpNotches = new List<int>();
                            List<IBioPolymerWithSetMods> allBmpPeptides = new List<IBioPolymerWithSetMods>();

                            foreach (var (Notch, Peptide) in psm.BestMatchingBioPolymersWithSetMods)
                            {
                                allBmpNotches.Add(Notch);
                                allBmpPeptides.Add(Peptide);
                                PsmData pd = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, Peptide, Notch, !Peptide.Parent.IsDecoy);
                                var pepValuePrediction = threadPredictionEngine.Predict(pd);
                                pepValuePredictions.Add(pepValuePrediction.Probability);
                                //A score is available using the variable pepvaluePrediction.Score
                            }

                            GetIndiciesOfPeptidesToRemove(indiciesOfPeptidesToRemove, pepValuePredictions);
                            int peptidesRemoved = 0;
                            RemoveBestMatchingPeptidesWithLowPEP(psm, indiciesOfPeptidesToRemove, allBmpNotches, allBmpPeptides, pepValuePredictions, ref peptidesRemoved);
                            ambigousPeptidesRemovedinThread += peptidesRemoved;

                            psm.PsmFdrInfo.PEP = 1 - pepValuePredictions.Max();
                            psm.PeptideFdrInfo.PEP = 1 - pepValuePredictions.Max();
                        }
                        
                    }

                    lock (lockObject)
                    {
                        ambiguousPeptidesResolved += ambigousPeptidesRemovedinThread;
                    }
                });

            return ambiguousPeptidesResolved;
        }

        public static IEnumerable<PsmData> CreatePsmData(string searchType,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            List<PeptideMatchGroup> allPeptides,
            List<int> peptideIndices,
            Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified,
            Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified,
            Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode)
        {
            object psmDataListLock = new object();
            List<PsmData> psmDataList = new List<PsmData>();
            List<double> psmOrder = new List<double>();
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, peptideIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<PsmData> localPsmDataList = new List<PsmData>();
                    List<double> localPsmOrder = new List<double>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        foreach (SpectralMatch psm in allPeptides[peptideIndices[i]])
                        {
                            // Stop loop if canceled
                            if (GlobalVariables.StopLoops) { return; }

                            PsmData newPsmData = new PsmData();
                            if (searchType == "crosslink")
                            {
                                CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)allPeptides[peptideIndices[i]].BestMatch;

                                bool label;
                                if (csm.IsDecoy || csm.BetaPeptide.IsDecoy)
                                {
                                    label = false;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                                }
                                else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && psm.FdrInfo.QValue <= FdrThresholdForPositiveExamples)
                                {
                                    label = true;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                                }
                                else
                                {
                                    continue;
                                }
                                localPsmDataList.Add(newPsmData);
                                localPsmOrder.Add(i);
                            }
                            else
                            {
                                // I think sometimes things get duplicated here but idk what to do about it
                                foreach(var (notch, peptideWithSetMods) in psm.BestMatchingBioPolymersWithSetMods)
                                {
                                    bool label;
                                    if (peptideWithSetMods.Parent.IsDecoy)
                                    {
                                        label = false;
                                        newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                    }
                                    else if (!peptideWithSetMods.Parent.IsDecoy && psm.FdrInfo.QValue <= FdrThresholdForPositiveExamples)
                                    {
                                        label = true;
                                        newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                    }
                                    else
                                    {
                                        continue;
                                    }
                                    if (newPsmData != null)
                                    {
                                        localPsmDataList.Add(newPsmData);
                                        localPsmOrder.Add(i);
                                    }
                                }
                            }
                        }
                    }

                    lock (psmDataListLock)
                    {
                        psmDataList.AddRange(localPsmDataList);
                        psmOrder.AddRange(localPsmOrder);
                    }
                });
            PsmData[] pda = psmDataList.ToArray();
            double[] order = psmOrder.ToArray();

            Array.Sort(order, pda);//this sorts both arrays thru sorting the array in position one. The order array, keeps track of the positon in the original psms list and returns the PsmData array in that same order.

            return pda.AsEnumerable();
        }

        public static IEnumerable<PsmData> CreatePsmDataFromPEP_QValue(string searchType,
           List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
           List<PeptideMatchGroup> allPeptides,
           List<int> peptideIndices,
           Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified,
           Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified,
           Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode)
        {
            object psmDataListLock = new object();
            List<PsmData> psmDataList = new List<PsmData>();
            List<double> psmOrder = new List<double>();
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, peptideIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<PsmData> localPsmDataList = new List<PsmData>();
                    List<double> localPsmOrder = new List<double>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        foreach (SpectralMatch psm in allPeptides[peptideIndices[i]])
                        {

                            // Stop loop if canceled
                            if (GlobalVariables.StopLoops) { return; }

                            PsmData newPsmData = new PsmData();
                            if (searchType == "crosslink")
                            {
                                CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)allPeptides[peptideIndices[i]].BestMatch;

                                bool label;
                                if (csm.IsDecoy || csm.BetaPeptide.IsDecoy)
                                {
                                    label = false;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                                }
                                else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && psm.PeptideFdrInfo.PEP_QValue <= FdrThresholdForPositiveExamples)
                                {
                                    label = true;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                                }
                                else
                                {
                                    continue;
                                }
                                localPsmDataList.Add(newPsmData);
                                localPsmOrder.Add(i);
                            }
                            else
                            {
                                double bmp = 0;
                                foreach (var (notch, peptideWithSetMods) in psm.BestMatchingBioPolymersWithSetMods)
                                {
                                    bool label;
                                    double bmpc = psm.BestMatchingBioPolymersWithSetMods.Count();
                                    if (peptideWithSetMods.Parent.IsDecoy)
                                    {
                                        label = false;
                                        newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                    }
                                    else if (!peptideWithSetMods.Parent.IsDecoy && psm.FdrInfo.PEP_QValue <= FdrThresholdForPositiveExamples)
                                    {
                                        label = true;
                                        newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                    }
                                    else
                                    {
                                        continue;
                                    }
                                    localPsmDataList.Add(newPsmData);
                                    localPsmOrder.Add(i + (bmp / bmpc / 2.0));
                                    bmp += 1.0;
                                }
                            }
                        }
                    }
                    lock (psmDataListLock)
                    {
                        psmDataList.AddRange(localPsmDataList);
                        psmOrder.AddRange(localPsmOrder);
                    }
                });
            PsmData[] pda = psmDataList.ToArray();
            double[] order = psmOrder.ToArray();

            Array.Sort(order, pda);//this sorts both arrays thru sorting the array in position one. The order array, keeps track of the positon in the original psms list and returns the PsmData array in that same order.

            return pda.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataEntry(string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, SpectralMatch psm, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, IBioPolymerWithSetMods selectedPeptide, int notchToUse, bool label)
        {
            double normalizationFactor = selectedPeptide.BaseSequence.Length;
            float totalMatchingFragmentCount = 0;
            float intensity = 0;
            float chargeDifference = 0;
            float deltaScore = 0;
            int notch = 0;
            float ambiguity = 0;
            float modCount = 0;
            float absoluteFragmentMassError = 0;

            float missedCleavages = 0;
            float longestSeq = 0;
            float complementaryIonCount = 0;
            float hydrophobicityZscore = float.NaN;
            bool isVariantPeptide = false;

            //crosslink specific features
            float alphaIntensity = 0;
            float betaIntensity = 0;
            float longestFragmentIonSeries_Alpha = 0;
            float longestFragmentIonSeries_Beta = 0;
            float isDeadEnd = 0;
            float isLoop = 0;
            float isInter = 0;
            float isIntra = 0;
            float spectralAngle = 0;
            float hasSpectralAngle = 0;

            if (searchType != "crosslink")
            {
                if (searchType == "top-down")
                {
                    normalizationFactor /= 10.0;
                }
                totalMatchingFragmentCount = (float)(Math.Round(psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide].Count / normalizationFactor * 10, 0));
                intensity = (float)Math.Min(50, Math.Round((psm.Score - (int)psm.Score) / normalizationFactor * 100.0, 0));
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                deltaScore = (float)Math.Round(psm.DeltaScore / normalizationFactor * 10.0, 0);
                notch = notchToUse;
                modCount = Math.Min((float)selectedPeptide.AllModsOneIsNterminus.Keys.Count(), 10);
                if (psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide]?.Count() > 0)
                {
                    absoluteFragmentMassError = (float)Math.Min(100.0, Math.Round(10.0 * Math.Abs(GetAverageFragmentMassError(psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide]) - fileSpecificMedianFragmentMassErrors[Path.GetFileName(psm.FullFilePath)])));
                }

                ambiguity = Math.Min((float)(psm.BioPolymersWithSetModsToMatchingFragments.Keys.Count - 1), 10);
                longestSeq = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(psm.BioPolymersWithSetModsToMatchingFragments, selectedPeptide) / normalizationFactor * 10, 0);
                complementaryIonCount = (float)Math.Round(SpectralMatch.GetCountComplementaryIons(psm.BioPolymersWithSetModsToMatchingFragments, selectedPeptide) / normalizationFactor * 10, 0);
                isVariantPeptide = PeptideIsVariant(selectedPeptide);
                spectralAngle = (float)psm.SpectralAngle;

                if (PsmHasSpectralAngle(psm))
                {
                    hasSpectralAngle = 1;
                }

                if (psm.DigestionParams.Protease.Name != "top-down")
                {
                    missedCleavages = selectedPeptide.MissedCleavages;
                    bool fileIsCzeSeparationType = fileSpecificParameters.Any(p => Path.GetFileName(p.fileName) == Path.GetFileName(psm.FullFilePath) && p.fileSpecificParameters.SeparationType == "CZE");

                    if (!fileIsCzeSeparationType)
                    {
                        if (selectedPeptide.BaseSequence.Equals(selectedPeptide.FullSequence))
                        {
                            hydrophobicityZscore = (float)Math.Round(GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_unmodified) * 10.0, 0);
                        }
                        else
                        {
                            hydrophobicityZscore = (float)Math.Round(GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_modified) * 10.0, 0);
                        }
                    }
                    else
                    {
                        hydrophobicityZscore = (float)Math.Round(GetMobilityZScore(psm, selectedPeptide) * 10.0, 0);
                    }
                }
                //this is not for actual crosslinks but for the byproducts of crosslink loop links, deadends, etc.
                if (psm is CrosslinkSpectralMatch)
                {
                    CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                    isDeadEnd = Convert.ToSingle((csm.CrossType == PsmCrossType.DeadEnd) || (csm.CrossType == PsmCrossType.DeadEndH2O) || (csm.CrossType == PsmCrossType.DeadEndNH2) || (csm.CrossType == PsmCrossType.DeadEndTris));
                    isLoop = Convert.ToSingle(csm.CrossType == PsmCrossType.Loop);
                }
            }
            else
            {
                CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                PeptideWithSetModifications selectedAlphaPeptide = csm.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide as PeptideWithSetModifications).First();
                PeptideWithSetModifications selectedBetaPeptide = csm.BetaPeptide?.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide as PeptideWithSetModifications).First();

                float alphaNormalizationFactor = selectedAlphaPeptide.BaseSequence.Length;
                float betaNormalizationFactor = selectedBetaPeptide == null ? (float)0 : selectedBetaPeptide.BaseSequence.Length;
                float totalNormalizationFactor = alphaNormalizationFactor + betaNormalizationFactor;

                totalMatchingFragmentCount = (float)Math.Round(csm.XLTotalScore / totalNormalizationFactor * 10, 0);

                //Compute fragment mass error
                int alphaCount = 0;
                float alphaError = 0;
                if (csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide]?.Count > 0)
                {
                    alphaCount = csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide].Count;
                    alphaError = Math.Abs(GetAverageFragmentMassError(csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide]));
                }
                int betaCount = 0;
                float betaError = 0;
                if (csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide]?.Count > 0)
                {
                    betaCount = csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide].Count;
                    betaError = Math.Abs(GetAverageFragmentMassError(csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide]));
                }

                float averageError = 0;
                if ((alphaCount + betaCount) > 0)
                {
                    averageError = (alphaCount * alphaError + betaCount * betaError) / (alphaCount + betaCount);
                }

                absoluteFragmentMassError = (float)Math.Min(100, Math.Round(averageError - fileSpecificMedianFragmentMassErrors[Path.GetFileName(csm.FullFilePath)] * 10.0, 0));
                //End compute fragment mass error

                deltaScore = (float)Math.Round(csm.DeltaScore / totalNormalizationFactor * 10.0, 0);
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                alphaIntensity = (float)Math.Min(100, Math.Round((csm.Score - (int)csm.Score) / alphaNormalizationFactor * 100.0, 0));
                betaIntensity = csm.BetaPeptide == null ? (float)0 : (float)Math.Min(100.0, Math.Round((csm.BetaPeptide.Score - (int)csm.BetaPeptide.Score) / betaNormalizationFactor * 100.0, 0));
                longestFragmentIonSeries_Alpha = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(csm.BioPolymersWithSetModsToMatchingFragments, selectedAlphaPeptide) / alphaNormalizationFactor * 10.0, 0);
                longestFragmentIonSeries_Beta = selectedBetaPeptide == null ? (float)0 : SpectralMatch.GetLongestIonSeriesBidirectional(csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments, selectedBetaPeptide) / betaNormalizationFactor;
                longestFragmentIonSeries_Beta = (float)Math.Round(longestFragmentIonSeries_Beta * 10.0, 0);
                isInter = Convert.ToSingle(csm.CrossType == PsmCrossType.Inter);
                isIntra = Convert.ToSingle(csm.CrossType == PsmCrossType.Intra);
            }

            psm.PsmData_forPEPandPercolator = new PsmData
            {
                TotalMatchingFragmentCount = totalMatchingFragmentCount,
                Intensity = intensity,
                PrecursorChargeDiffToMode = chargeDifference,
                DeltaScore = deltaScore,
                Notch = notch,
                ModsCount = modCount,
                AbsoluteAverageFragmentMassErrorFromMedian = absoluteFragmentMassError,
                MissedCleavagesCount = missedCleavages,
                Ambiguity = ambiguity,
                LongestFragmentIonSeries = longestSeq,
                ComplementaryIonCount = complementaryIonCount,
                HydrophobicityZScore = hydrophobicityZscore,
                IsVariantPeptide = Convert.ToSingle(isVariantPeptide),

                AlphaIntensity = alphaIntensity,
                BetaIntensity = betaIntensity,
                LongestFragmentIonSeries_Alpha = longestFragmentIonSeries_Alpha,
                LongestFragmentIonSeries_Beta = longestFragmentIonSeries_Beta,
                IsDeadEnd = isDeadEnd,
                IsLoop = isLoop,
                IsInter = isInter,
                IsIntra = isIntra,

                Label = label,

                SpectralAngle = spectralAngle,
                HasSpectralAngle = hasSpectralAngle
            };

            return psm.PsmData_forPEPandPercolator;
        }

        public static List<int>[] GetPeptideGroupIndices(List<PeptideMatchGroup> peptides, int numGroups)
        {
            List<int>[] groupsOfIndices = new List<int>[numGroups];

            List<int> targetIndices = new List<int>();
            List<int> decoyIndices = new List<int>();
            for(int i = 0; i< peptides.Count; i++)
            {
                if (peptides[i].Any(p => p.IsDecoy))
                {
                   decoyIndices.Add(i);
                }
                else
                {
                   targetIndices.Add(i);
                }
            }

            var targetIndexGroups = DivideListIntoGroups(targetIndices, numGroups);
            var decoyIndexGroups = DivideListIntoGroups(decoyIndices, numGroups);

            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndices[i] = targetIndexGroups[i].Concat(decoyIndexGroups[i]).ToList();
            }

            return groupsOfIndices;
        }

        /// <summary>
        /// This takes in a list of ints, and partitions them into numGroups partitions,
        /// e.g., partition 1 = [0, 4, 8...], partition 2 = [1, 5, 9...], etc.
        /// </summary>
        /// <returns>A list containing numGroups partitions (lists of ints) </returns>
        static List<List<int>> DivideListIntoGroups(List<int> list, int numGroups)
        {
            var groups = new List<List<int>>();
            for (int i = 0; i < numGroups; i++)
            {
                groups.Add(new List<int>());
            }

            int mainIndex = 0;
            while (mainIndex < list.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && mainIndex < list.Count)
                {
                    groups[subIndex].Add(mainIndex);

                    subIndex++;
                    mainIndex++;
                }
            }

            return groups;
        }

        public static void RemoveBestMatchingPeptidesWithLowPEP(SpectralMatch psm, List<int> indiciesOfPeptidesToRemove, List<int> notches, List<IBioPolymerWithSetMods> pwsmList, List<double> pepValuePredictions, ref int ambiguousPeptidesRemovedCount)
        {
            foreach (int i in indiciesOfPeptidesToRemove)
            {
                psm.RemoveThisAmbiguousPeptide(notches[i], pwsmList[i]);
                ambiguousPeptidesRemovedCount++;
            }
            psm.PsmFdrInfo.PEP = 1 - pepValuePredictions.Max();
            psm.PeptideFdrInfo.PEP = 1 - pepValuePredictions.Max();
        }

        /// <summary>
        /// Given a set of PEP values, this method will find the indicies of BestMatchingBioPolymersWithSetMods that are not within the required tolerance
        /// This method will also remove the low scoring predictions from the set.
        /// </summary>
        public static void GetIndiciesOfPeptidesToRemove(List<int> indiciesOfPeptidesToRemove, List<double> pepValuePredictions)
        {
            double highestPredictedPEPValue = pepValuePredictions.Max();
            for (int i = 0; i < pepValuePredictions.Count; i++)
            {
                if ((highestPredictedPEPValue - pepValuePredictions[i]) > AbsoluteProbabilityThatDistinguishesPeptides)
                {
                    indiciesOfPeptidesToRemove.Add(i);
                }
            }

            foreach (int i in indiciesOfPeptidesToRemove.OrderByDescending(p => p))
            {
                pepValuePredictions.RemoveAt(i);
            }
        }


        #region Dictionary Construction, scoring helpers
        /// <summary>
        /// Here we're getting the most common charge state for precursors that are Targets with q<=0.01.

        public static int GetChargeStateMode(IEnumerable<SpectralMatch> psms)
        {
            return psms.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.01).Select(p => p.ScanPrecursorCharge).GroupBy(n => n).OrderByDescending(g => g.Count()).Select(g => g.Key).FirstOrDefault();
        }

        public static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeHydrophobicityValues(IEnumerable<SpectralMatch> psms, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, bool computeHydrophobicitiesforModifiedPeptides)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            //TODO change the tuple so the values have names
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtHydrophobicityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = fileSpecificParameters.Where(s => s.fileSpecificParameters.SeparationType == "HPLC").Select(s => Path.GetFileName(s.fileName)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrophobicities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach ((int notch, IBioPolymerWithSetMods pwsm) in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(pwsm.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(pwsm.FullSequence);

                        double predictedHydrophobicity = pwsm is PeptideWithSetModifications pep ? calc.ScoreSequence(pep) : 0;

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        //First block of if statement is for modified peptides.
                        if (pwsm.AllModsOneIsNterminus.Any() && computeHydrophobicitiesforModifiedPeptides)
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                        //this second block of if statment is for unmodified peptides.
                        else if (!pwsm.AllModsOneIsNterminus.Any() && !computeHydrophobicitiesforModifiedPeptides)
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                    }
                }

                List<double> allSquaredHyrophobicityDifferences = new List<double>();

                foreach (int retentionTimeBin in hydrophobicities.Keys)
                {
                    //TODO consider using inner-quartile range instead of standard deviation
                    double averageHydrophobicity = hydrophobicities[retentionTimeBin].Average();
                    averagesCommaStandardDeviations.Add(retentionTimeBin, new Tuple<double, double>(averageHydrophobicity, hydrophobicities[retentionTimeBin].StandardDeviation()));
                    foreach (double hydrophobicity in hydrophobicities[retentionTimeBin])
                    {
                        double difference = Math.Abs(hydrophobicity - averageHydrophobicity);
                        if (!double.IsNaN(difference) && difference > 0)
                        {
                            allSquaredHyrophobicityDifferences.Add(Math.Pow(difference, 2));
                        }
                    }
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.
                double globalStDev = 1;
                if (allSquaredHyrophobicityDifferences.Count() > 1)
                {
                    globalStDev = Math.Sqrt(allSquaredHyrophobicityDifferences.Sum() / (allSquaredHyrophobicityDifferences.Count() - 1));
                }

                Dictionary<int, Tuple<double, double>> stDevsToChange = new Dictionary<int, Tuple<double, double>>();
                foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
                {
                    //add stability. not allowing stdevs that are too small or too large at one position relative to the global stdev
                    //here we are finding which stdevs are out of whack.
                    if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.5 || (item.Value.Item2 / globalStDev) > 3)
                    {
                        Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                        stDevsToChange.Add(item.Key, pair);
                    }
                }
                //here we are replacing the stdevs that are out of whack.
                foreach (int key in stDevsToChange.Keys)
                {
                    averagesCommaStandardDeviations[key] = stDevsToChange[key];
                }

                rtHydrophobicityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtHydrophobicityAvgDev;
        }

        public static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeMobilityValues(IEnumerable<SpectralMatch> psms, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters)
        {
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtMobilityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = fileSpecificParameters.Where(s => s.fileSpecificParameters.SeparationType == "CZE").Select(s => Path.GetFileName(s.fileName)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> mobilities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach ((int notch, IBioPolymerWithSetMods pwsm) in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(pwsm.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(pwsm.FullSequence);

                        double predictedMobility = pwsm is PeptideWithSetModifications pep ? 100.0 * GetCifuentesMobility(pep) : 0;

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        if (mobilities.ContainsKey(possibleKey))
                        {
                            mobilities[possibleKey].Add(predictedMobility);
                        }
                        else
                        {
                            mobilities.Add(possibleKey, new List<double> { predictedMobility });
                        }
                    }
                }

                List<double> allSquaredMobilityDifferences = new List<double>();

                foreach (int retentionTimeBin in mobilities.Keys)
                {
                    //TODO consider using inner-quartile range instead of standard deviation
                    double averageMobility = mobilities[retentionTimeBin].Average();
                    averagesCommaStandardDeviations.Add(retentionTimeBin, new Tuple<double, double>(averageMobility, mobilities[retentionTimeBin].StandardDeviation()));
                    foreach (double hydrophobicity in mobilities[retentionTimeBin])
                    {
                        double difference = Math.Abs(hydrophobicity - averageMobility);
                        if (!double.IsNaN(difference) && difference > 0)
                        {
                            allSquaredMobilityDifferences.Add(Math.Pow(difference, 2));
                        }
                    }
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.
                double globalStDev = 1;
                if (allSquaredMobilityDifferences.Count() > 1)
                {
                    globalStDev = Math.Sqrt(allSquaredMobilityDifferences.Sum() / (allSquaredMobilityDifferences.Count() - 1));
                }

                Dictionary<int, Tuple<double, double>> stDevsToChange = new Dictionary<int, Tuple<double, double>>();

                GetStDevsToChange(stDevsToChange, averagesCommaStandardDeviations, globalStDev);
                UpdateOutOfRangeStDevsWithGlobalAverage(stDevsToChange, averagesCommaStandardDeviations);

                rtMobilityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtMobilityAvgDev;
        }

        /// <summary>
        /// This gathers a set of standard deviations that are outside the range of acceptable.
        /// </summary>
        public static void GetStDevsToChange(Dictionary<int, Tuple<double, double>> stDevsToChange, Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations, double globalStDev)
        {
            foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
            {
                //add stability. not allowing stdevs that are too small or too large at one position relative to the global stdev
                //here we are finding which stdevs are out of whack.
                if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.05 || (item.Value.Item2 / globalStDev) > 3)
                {
                    Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                    stDevsToChange.Add(item.Key, pair);
                }
            }
        }

        /// <summary>
        /// here we are replacing the stdevs that are out of whack.
        /// </summary>
        public static void UpdateOutOfRangeStDevsWithGlobalAverage(Dictionary<int, Tuple<double, double>> stDevsToChange, Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations)
        {
            foreach (int key in stDevsToChange.Keys)
            {
                averagesCommaStandardDeviations[key] = stDevsToChange[key];
            }
        }

        private static double GetCifuentesMobility(IBioPolymerWithSetMods pwsm)
        {
            int charge = 1 + pwsm.BaseSequence.Count(f => f == 'K') + pwsm.BaseSequence.Count(f => f == 'R') + pwsm.BaseSequence.Count(f => f == 'H') - CountModificationsThatShiftMobility(pwsm.AllModsOneIsNterminus.Values.AsEnumerable());// the 1 + is for N-terminal

            double mobility = (Math.Log(1 + 0.35 * (double)charge)) / Math.Pow(pwsm.MonoisotopicMass, 0.411);

            return mobility;
        }

        private static float GetSSRCalcHydrophobicityZScore(SpectralMatch psm, IBioPolymerWithSetMods Peptide, Dictionary<string, Dictionary<int, Tuple<double, double>>> d)
        {
            //Using SSRCalc3 but probably any number of different calculators could be used instead. One could also use the CE mobility.
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            double hydrophobicityZscore = double.NaN;

            if (d.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (d[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedHydrophobicity = Peptide is PeptideWithSetModifications pep ? calc.ScoreSequence(pep) : 0;

                    hydrophobicityZscore = Math.Abs(d[Path.GetFileName(psm.FullFilePath)][time].Item1 - predictedHydrophobicity) / d[Path.GetFileName(psm.FullFilePath)][time].Item2;
                }
            }

            double maxHydrophobicityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
            if (double.IsNaN(hydrophobicityZscore) || double.IsInfinity(hydrophobicityZscore) || hydrophobicityZscore > maxHydrophobicityZscore)
            {
                hydrophobicityZscore = maxHydrophobicityZscore;
            }

            return (float)hydrophobicityZscore;
        }

        private static float GetMobilityZScore(SpectralMatch psm, IBioPolymerWithSetMods selectedPeptide)
        {
            double mobilityZScore = double.NaN;

            if (fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedMobility = 100.0 * GetCifuentesMobility(selectedPeptide);

                    mobilityZScore = Math.Abs(fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item1 - predictedMobility) / fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item2;
                }
            }

            double maxMobilityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
            if (double.IsNaN(mobilityZScore) || double.IsInfinity(mobilityZScore) || mobilityZScore > maxMobilityZscore)
            {
                mobilityZScore = maxMobilityZscore;
            }

            return (float)mobilityZScore;
        }

        private static bool PeptideIsVariant(IBioPolymerWithSetMods bpwsm)
        {
            if (bpwsm is not PeptideWithSetModifications pwsm)
                return false;

            bool identifiedVariant = false;
            if (pwsm.Protein.AppliedSequenceVariations.Count() > 0)
            {
                foreach (var variant in pwsm.Protein.AppliedSequenceVariations)
                {
                    if (pwsm.IntersectsAndIdentifiesVariation(variant).identifies)
                    {
                        identifiedVariant = true;
                        break;
                    }
                }
            }
            return identifiedVariant;
        }

        private static bool PsmHasSpectralAngle(SpectralMatch psm)
        {
            return psm.SpectralAngle >= 0;
        }

        public static bool ContainsModificationsThatShiftMobility(IEnumerable<Modification> modifications)
        {
            List<string> shiftingModifications = new List<string> { "Acetylation", "Ammonia loss", "Carbamyl", "Deamidation", "Formylation",
                "N2-acetylarginine", "N6-acetyllysine", "N-acetylalanine", "N-acetylaspartate", "N-acetylcysteine", "N-acetylglutamate", "N-acetylglycine",
                "N-acetylisoleucine", "N-acetylmethionine", "N-acetylproline", "N-acetylserine", "N-acetylthreonine", "N-acetyltyrosine", "N-acetylvaline",
                "Phosphorylation", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine", "Sulfonation" };

            return shiftingModifications.Concat(modifications.Select(m => m.OriginalId).Distinct()).GroupBy(s => s).Where(s => s.Count() > 1).Any();
        }

        public static int CountModificationsThatShiftMobility(IEnumerable<Modification> modifications)
        {
            List<string> shiftingModifications = new List<string> { "Acetylation", "Ammonia loss", "Carbamyl", "Deamidation", "Formylation",
                "N2-acetylarginine", "N6-acetyllysine", "N-acetylalanine", "N-acetylaspartate", "N-acetylcysteine", "N-acetylglutamate", "N-acetylglycine",
                "N-acetylisoleucine", "N-acetylmethionine", "N-acetylproline", "N-acetylserine", "N-acetylthreonine", "N-acetyltyrosine", "N-acetylvaline",
                "Phosphorylation", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine", "Sulfonation" };

            return modifications.Select(n => n.OriginalId).Intersect(shiftingModifications).Count();
        }

        public static Dictionary<string, float> GetFileSpecificMedianFragmentMassError(IEnumerable<SpectralMatch> psms)
        {
            Dictionary<string, float> fileSpecificMassErrors = new Dictionary<string, float>();
            foreach (string filename in psms.Select(p => Path.GetFileName(p.FullFilePath)).Distinct())
            {
                fileSpecificMassErrors.Add(filename, GetMedianAverageMassError(psms.Where(p => Path.GetFileName(p.FullFilePath) == filename)));
            }
            return fileSpecificMassErrors;
        }

        public static float GetMedianAverageMassError(IEnumerable<SpectralMatch> psms)
        {
            List<float> averageMassErrors = new List<float>();
            foreach (SpectralMatch psm in psms)
            {
                {
                    foreach (KeyValuePair<IBioPolymerWithSetMods, List<MatchedFragmentIon>> peptide_MFI in psm.BioPolymersWithSetModsToMatchingFragments)
                    {
                        if (peptide_MFI.Value != null && peptide_MFI.Value.Count > 0)
                        {
                            averageMassErrors.Add(GetAverageFragmentMassError(peptide_MFI.Value));
                        }
                    }
                }
            }
            return averageMassErrors.Median();
        }

        public static float GetAverageFragmentMassError(IEnumerable<MatchedFragmentIon> matchedIons)
        {
            var matchedIonsGroupedByProductType = matchedIons.GroupBy(i => i.NeutralTheoreticalProduct.ProductType).OrderBy(i => i.Key).ToList();
            List<float> massErrors = new List<float>();
            foreach (var productType in matchedIonsGroupedByProductType)
            {
                var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.FragmentNumber)
                    .ToList();

                for (int i = 0; i < products.Count; i++)
                {
                    MatchedFragmentIon ion = products[i];

                    float massError = (float)(ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass);
                    float ppmMassError = (float)(massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6);
                    massErrors.Add(ppmMassError);
                }
            }

            return massErrors.Average();
        }

        #endregion
    }
}
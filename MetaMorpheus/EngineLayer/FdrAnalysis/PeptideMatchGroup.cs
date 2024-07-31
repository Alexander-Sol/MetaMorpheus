﻿using Omics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class PeptideMatchGroup : IEnumerable<SpectralMatch>
    {
        public string PeptideFullSequence { get; }
        public List<SpectralMatch> SpectralMatches { get; }

        /// <summary>
        /// This class groups all spectra associated with a given peptide together,
        /// to facilitate the calculation of PEP values.
        /// </summary>
        /// <param name="fullPeptideSeq"> The full sequence to be used for grouping</param>
        /// <param name="spectralMatches"> Every spectral match that matches the full sequence</param>
        public PeptideMatchGroup(string fullPeptideSeq, List<SpectralMatch> spectralMatches)
        {
            PeptideFullSequence = fullPeptideSeq;
            SpectralMatches = spectralMatches;
        }

        public static List<PeptideMatchGroup> GroupByFullSequence(List<SpectralMatch> spectralMatches)
        {
            // This groups psms by full sequences. All ambiguous PSMs are grouped together (full sequence == null), this prevents accidental cross contamination during training.
            return spectralMatches.GroupBy(p => p.FullSequence)
                .Select(group => new PeptideMatchGroup(group.Key, group.ToList()))
                .OrderByDescending(matchGroup => matchGroup.Count())
                .ThenByDescending(matchGroup => matchGroup.BestMatch.Score)
                .ToList();
        }

        public static List<PeptideMatchGroup> GroupByBaseSequence(List<SpectralMatch> spectralMatches)
        {
            // This groups psms by base sequence, ensuring that PSMs with the same base sequence but different modifications are grouped together when training.
            return spectralMatches.GroupBy(p => p.BaseSequence)
                .Select(group => new PeptideMatchGroup(group.Key, group.ToList()))
                .OrderByDescending(matchGroup => matchGroup.Count())
                .ThenByDescending(matchGroup => matchGroup.BestMatch.Score)
                .ToList();
        }

        public IEnumerable<SpectralMatch> GetBestMatchByMod()
        {
            return SpectralMatches.GroupBy(p => p.FullSequence).Select(g => g.MaxBy(p => p));
        }

        /// <summary>
        /// This function is called if there aren't enough peptides to train at the peptide level
        /// </summary>
        /// <param name="spectralMatches"></param>
        /// <returns></returns>
        public static List<PeptideMatchGroup> GroupByIndividualPsm(List<SpectralMatch> spectralMatches)
        {
            return spectralMatches.Select(psm => new PeptideMatchGroup(psm.FullSequence, new List<SpectralMatch> { psm }))
                .ToList();
        }

        public SpectralMatch BestMatch => SpectralMatches.MaxBy(match => match);

        public IEnumerator<SpectralMatch> GetEnumerator()
        {
            return SpectralMatches.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

    }
}
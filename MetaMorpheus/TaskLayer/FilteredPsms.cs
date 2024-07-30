﻿using Easy.Common.Extensions;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public enum FilterType
    {
        QValue,
        PepQValue
    }

    /// <summary>
    /// Contains a filtered list of PSMs
    /// </summary>
    public class FilteredPsms : IEnumerable<SpectralMatch>
    {
        public List<SpectralMatch> Psms { get; set; }
        /// <summary>
        /// Filter type can have only two values: "q-value" or "pep q-value"
        /// </summary>
        public FilterType FilterType { get; }
        public double FilterThreshold { get; }
        public bool FilteringNotPerformed { get; }
        public bool PeptideLevelFiltering { get; }
        public FilteredPsms(List<SpectralMatch> psms, FilterType filterType, double filterThreshold, bool filteringNotPerformed, bool peptideLevelFiltering)
        {
            Psms = psms;
            FilterType = filterType;
            FilterThreshold = filterThreshold;
            FilteringNotPerformed = filteringNotPerformed;
            PeptideLevelFiltering = peptideLevelFiltering;
        }

        private bool AboveThreshold(SpectralMatch psm)
        {
            if (psm.GetFdrInfo(PeptideLevelFiltering) == null) return false;

            switch (FilterType)
            {
                case FilterType.PepQValue:
                    return psm.GetFdrInfo(PeptideLevelFiltering).PEP_QValue <= FilterThreshold;
                default:
                    return psm.GetFdrInfo(PeptideLevelFiltering).QValue <= FilterThreshold && psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch <= FilterThreshold;
            }
        }

        public string GetFilterTypeString()
        {
            return FilterType == FilterType.PepQValue ? "pep q-value" : "q-value";
        }

        /// <summary>
        /// Returns the number of PSMs that passed the filtering criteria
        /// </summary>
        public int PsmsAboveThreshold => Psms.Count(psm => AboveThreshold(psm));

        public IEnumerator<SpectralMatch> GetEnumerator()
        {
            return Psms.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return Psms.GetEnumerator();
        }
    }
}

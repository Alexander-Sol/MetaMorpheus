using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.IO;
using UsefulProteomicsDatabases;

namespace EngineLayer.SequenceCoverage
{
    internal class ProteaseCoverage
    {
        public readonly Protease Protease;
        public readonly Protein Protein;
        public readonly Char[] Residues;
        public bool[] PeptideCoverage { get; private set; }
        public int[] FragmentCoverage { get; private set; }

        public ProteaseCoverage()
        {

        }
    }
}

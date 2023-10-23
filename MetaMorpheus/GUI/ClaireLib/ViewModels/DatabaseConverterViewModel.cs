﻿using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
using EngineLayer;
using GuiFunctions;
using pepXML.Generated;
using Proteomics;
using Readers;
using TopDownProteomics.IO.MzIdentMl;
using UsefulProteomicsDatabases;

namespace MetaMorpheusGUI
{
    public class DatabaseConverterViewModel : BaseViewModel
    {
        public ObservableCollection<string> DatabasePaths { get; }
        private string selectedDatbase;
        public string SelectedDatabase
        {
            get => selectedDatbase;
            set { selectedDatbase = value; OnPropertyChanged(nameof(SelectedDatabase)); }
        }

        public ObservableCollection<string> SearchResultPaths { get; }
        private string selectedSearchResult;
        public string SelectedSearchResult
        {
            get => selectedSearchResult;
            set { selectedSearchResult = value; OnPropertyChanged(nameof(SelectedSearchResult)); }
        }


        private string _outputDatabasePath;
        public string OutputDatabasePath
        {
            get => _outputDatabasePath ??= DatabasePaths.FirstOrDefault()?.Replace(".xml", ".fasta") ?? null;
            set { _outputDatabasePath = value; OnPropertyChanged(nameof(OutputDatabasePath)); }
        }

        #region Parameters

        private bool _generateDecoys;
        public bool GenerateDecoys
        {
            get => _generateDecoys;
            set
            {
                _generateDecoys = value;
                if (!value)
                    SelectedDecoyType = DecoyType.None;
                OnPropertyChanged(nameof(GenerateDecoys));
            }
        }

        private bool _generateTargets;
        public bool GenerateTargets
        {
            get => _generateTargets;
            set { _generateTargets = value; OnPropertyChanged(nameof(GenerateTargets)); }
        }

        public DecoyType[] DecoyTypes { get; set; }

        private DecoyType selectedDecoyType;
        public DecoyType SelectedDecoyType
        {
            get => selectedDecoyType;
            set { selectedDecoyType = value; OnPropertyChanged(nameof(SelectedDecoyType)); }
        }

        private bool _appendSearchResults;
        public bool AppendSearchResults
        {
            get => _appendSearchResults;
            set { _appendSearchResults = value; OnPropertyChanged(nameof(AppendSearchResults)); }
        }

        

        private bool _filterToFdr;
        public bool FilterToFdr
        {
            get => _filterToFdr;
            set { _filterToFdr = value; OnPropertyChanged(nameof(FilterToFdr)); }
        }

        private double _fdrCutoff;
        public double FdrCutoff
        {
            get => _fdrCutoff;
            set { _fdrCutoff = value; OnPropertyChanged(nameof(FdrCutoff)); }
        }

        //private bool


        // not implemented, only fasta can be exported
        public string[] OutputTypes { get; set; }
        private string selectedOutputType;
        public string SelectedOutputType
        {
            get => selectedOutputType;
            set { selectedOutputType = value; OnPropertyChanged(nameof(SelectedOutputType)); }
        }

        #endregion



        public DatabaseConverterViewModel()
        {
            DatabasePaths = new ObservableCollection<string>();
            SearchResultPaths = new ObservableCollection<string>();
            DecoyTypes = Enum.GetValues<DecoyType>()
                .Where(p => p != DecoyType.Random)
                .ToArray();
            SelectedDecoyType = DecoyType.None;
            OutputTypes = new[] { "fasta", "xml" };
            SelectedOutputType = "fasta";
            GenerateTargets = true;
            GenerateDecoys = false;
            AppendSearchResults = true;
            FilterToFdr = true;
            FdrCutoff = 0.01;

            RemoveDatabaseCommand = new RelayCommand(RemoveDatabaseFromDatabasePaths);
            RemoveSearchResultCommand = new RelayCommand(RemoveSearchResultFromSearchResultPaths);
            CreateSingleDatabaseCommand = new RelayCommand(CreateSingleDatabase);
            ClearDataCommand = new RelayCommand(Clear);
        }

        public ICommand RemoveDatabaseCommand { get; set; }
        private void RemoveDatabaseFromDatabasePaths()
        {
            DatabasePaths.Remove(SelectedDatabase);
            SelectedDatabase = DatabasePaths.FirstOrDefault();
        }

        public ICommand RemoveSearchResultCommand { get; set; }

        private void RemoveSearchResultFromSearchResultPaths()
        {
            SearchResultPaths.Remove(SelectedSearchResult);
            SelectedSearchResult = SearchResultPaths.FirstOrDefault();
        }

        public ICommand ClearDataCommand { get; set; }
        public void Clear()
        {
            DatabasePaths.Clear();
            SearchResultPaths.Clear();
            OutputDatabasePath = null;
        }

        private string GetFinalPath(string path)
        {
            // check if a file with this name already exists, if so add a number to the end within parenthesis. If that file still exists, increment the number by one and try again
            int fileCount = 1;
            string finalPath = path;
            while (System.IO.File.Exists(finalPath))
            {
                finalPath = path.Replace(".fasta", $"({fileCount}).fasta");
                fileCount++;
            }
            return finalPath;
        }

        public ICommand CreateSingleDatabaseCommand { get; set; }

        private void CreateSingleDatabase()
        {
            List<string> readingErrors = new();

            // load proteins from input database
            List<Protein> proteins = new List<Protein>();
            foreach (var database in DatabasePaths)
            {
                try
                {
                    if (database.EndsWith(".xml"))
                    {
                        proteins.AddRange(ProteinDbLoader.LoadProteinXML(database, GenerateTargets, SelectedDecoyType,
                            GlobalVariables.AllModsKnown, false,
                            new List<string>(), out Dictionary<string, Modification> unknownMods));
                        readingErrors.AddRange(unknownMods.Select(keyValuePair =>
                            $"unknown modificaiton found on {keyValuePair.Key} with type {keyValuePair.Value.IdWithMotif}"));
                    }
                    else if (database.EndsWith(".fasta"))
                    {
                        proteins.AddRange(ProteinDbLoader.LoadProteinFasta(database, GenerateTargets, SelectedDecoyType,
                            false, out List<string> error));
                        readingErrors.AddRange(error);
                    }
                }
                catch (Exception e)
                {
                    MessageBox.Show($"Error Reading in Database {database}\n{e.Message}");
                }
            }

            // write proteins to output database
            string finalPath = GetFinalPath(OutputDatabasePath);
            WriteDatabase(proteins, finalPath);

            // load in search results and collect fasta headers
            List<string> fastaLinesToAdd = new List<string>();
            foreach (var resultPath in SearchResultPaths)
            {
                try
                {
                    if (resultPath.EndsWith(".psmtsv"))
                    {
                        var psms = PsmTsvReader.ReadTsv(resultPath, out List<string> warnings);
                        readingErrors.AddRange(warnings);

                        if (FilterToFdr)
                            psms = psms.Where(p => p.QValue <= FdrCutoff).ToList();

                        foreach (var psm in psms)
                        {
                            fastaLinesToAdd.Add(psm.GetUniprotHeaderFromPsmFromTsv());
                            fastaLinesToAdd.Add(psm.BaseSeq);
                        }
                    }
                    else if (resultPath.ParseFileType().ToString().Contains("Toppic"))
                    {
                        ToppicSearchResultFile file = new ToppicSearchResultFile(resultPath);
                        var psms = file.Results;

                        if (FilterToFdr)
                            psms = psms.Where(p => p.EValue <= FdrCutoff).ToList();

                        foreach (var psm in psms)
                        {
                            fastaLinesToAdd.Add(psm.GetUniprotHeaderFromToppicPrsm());
                            fastaLinesToAdd.Add(psm.BaseSequence);
                        }


                         }
                }
                catch (Exception e)
                {
                    MessageBox.Show($"Error Reading in Search Result {resultPath}\n{e.Message}");
                }
            }

            // append fasta lines from search results to output database

            
            MessageBox.Show($"New Database Outputted to {finalPath}");
        }


        private void WriteDatabase(List<Protein> inputDatabaseProteins, string finalpath)
        {
            try
            {
                switch (SelectedOutputType)
                {
                    case "fasta":
                        ProteinDbWriter.WriteFastaDatabase(inputDatabaseProteins, finalpath, ">");
                        break;
                    case "xml":
                        var modDict = CreateModDictionary(inputDatabaseProteins);
                        ProteinDbWriter.WriteXmlDatabase(modDict, inputDatabaseProteins, finalpath);
                        break;
                    default:
                        throw new ArgumentException("not a valid database output type");
                }
            }
            catch (Exception e)
            {
                MessageBox.Show($"Database Writing Error: {e.Message}");
            }
        }

        private void AppendSearchResultsToDatabase(string finalPath, List<string> linesToAppend)
        {
            using var sw = new StreamWriter(finalPath, true);
            linesToAppend.ForEach(p => sw.WriteLine(p));
        }





        private Dictionary<string, HashSet<Tuple<int, Modification>>> CreateModDictionary(List<Protein> proteins)
        {
            Dictionary<string, HashSet<Tuple<int, Modification>>> modDict = new();
            HashSet<Tuple<int, Modification>> mods = new();
            foreach (var protein in proteins)
            {
                mods.Clear();
                foreach (var modificationDictEntry in protein.OneBasedPossibleLocalizedModifications)
                {
                    foreach (var mod in modificationDictEntry.Value)
                    {
                        var tuple = new Tuple<int, Modification>(modificationDictEntry.Key, mod);
                        if (modDict.TryGetValue(protein.Accession, out var hash))
                            hash.Add(tuple);
                        
                        else
                            modDict[protein.Accession] = new HashSet<Tuple<int, Modification>>() { tuple };
                    }
                }
            }

            return modDict;
        }

        internal void FileDropped(string path)
        {
            // database
            if (path.EndsWith(".xml") || path.EndsWith(".fasta"))
            {
                if (DatabasePaths.Contains(path))
                    return;

                DatabasePaths.Add(path);
                OnPropertyChanged(nameof(OutputDatabasePath));
            }
            else if (path.EndsWith(".psmtsv") || path.EndsWith(".tsv"))
            {
                if (SearchResultPaths.Contains(path))
                    return;

                SearchResultPaths.Add(path);
            }
            else
            {
                MessageBox.Show($"Cannot determine file type of and ignored file: {path}");
            }
        }
        
    }

    public static class IdExtensions
    {
        public static string GetUniprotHeaderFromPsmFromTsv(this PsmFromTsv psm)
        {
            var gene = "";
            var geneName = psm.GeneName.Split(',');
            if (geneName.Any())
                gene = geneName.First().Split(':')[1];

            
            var str =  $"mz|{psm.ProteinAccession}|{psm.ProteinName} OS={psm.OrganismName} GN={gene}";
            return str;
        }

        public static string GetUniprotHeaderFromToppicPrsm(this ToppicPrsm prms)
        {

        }
    }
}

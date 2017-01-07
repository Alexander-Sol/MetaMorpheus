﻿using MetaMorpheus;
using System.Collections.ObjectModel;
using System.Windows;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        private ObservableCollection<ModList> modFileList;
        private MyCalibrateTask myCalibrateTask;

        public CalibrateTaskWindow()
        {
            InitializeComponent();
        }

        public CalibrateTaskWindow(ObservableCollection<ModList> modFileList)
        {
            this.modFileList = modFileList;
        }

        public CalibrateTaskWindow(MyCalibrateTask myCalibrateTask, ObservableCollection<ModList> modFileList)
        {
            this.myCalibrateTask = myCalibrateTask;
            this.modFileList = modFileList;
        }

        internal MyTask TheTask { get; set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
        }
    }
}
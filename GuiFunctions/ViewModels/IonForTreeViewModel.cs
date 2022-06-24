﻿using GuiFunctions;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GuiFunctions
{
    public class IonForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        protected string _selectedColor;
        protected SolidColorBrush _colorBrush;

        #endregion

        #region Public Properties

        public ProductType IonType { get; set; }    
        public string IonName { get; set; }
        public string SelectedColor
        {
            get { return _selectedColor; }
            set 
            { 
                _selectedColor = value;
                ColorBrush = DrawnSequence.ParseColorBrushFromName(_selectedColor);
                OnPropertyChanged(nameof(SelectedColor));
            }
        }
        public SolidColorBrush ColorBrush
        {
            get { return _colorBrush; }
            set
            {
                _colorBrush = value;
                OnPropertyChanged(nameof(ColorBrush));
            }
        }

        public bool IsBeta { get; set; }
        public bool HasChanged { get; set; } = false;

        #endregion

        #region Constructor

        public IonForTreeViewModel(ProductType type, bool beta)
        {
            IonType = type;
            IonName = IonType.ToString() + " - Ion";
            IsBeta = beta;
            OxyColor color;
            if (IsBeta)
                color = MetaDrawSettings.BetaProductTypeToColor[IonType];
            else
                color = MetaDrawSettings.ProductTypeToColor[IonType];
            SelectedColor = AddSpaces(color.GetColorName());
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }

    }
}

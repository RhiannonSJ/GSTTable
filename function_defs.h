/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : June 2017
 *
 *
*/
#ifndef FUNCTION_DEFS
#define FUNCTION_DEFS

#include <fstream>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TObjArray.h"

using namespace std;

// -------------------------------------------------------------------------
//                          Typedefs 
// -------------------------------------------------------------------------

typedef map< string, vector < double > > m_outer; 
typedef map< string, double > m_map;

// -------------------------------------------------------------------------
// Normalisation:
// 
// Scale the event rate to the relevant nominal POT for the detector
//
// The arguments passed to the normalisation function are:
//      nEvents : The number of events we have generated
//      nPOT    : The corresponding POT, accessed from gst->GetWeight()
//      det     : Detector enumeration for nominal POT
//                0 = SBND
//                1 = MicroBooNE
//                2 = ICARUS
//
// -------------------------------------------------------------------------
double Norm ( const int nEvents,
              const double nPOT,
              const int det);

// -------------------------------------------------------------------------
// Make a map of the number of final state particles in each model
// configuration
// -------------------------------------------------------------------------
void FSPNumbers( TTree *event_tree,
                 m_map &n_fsp );

// -------------------------------------------------------------------------
// Make a map of the number of different final state interactions 
// in each model configuration
// -------------------------------------------------------------------------
void FSINumbers( TTree *event_tree,
                 double norm,
                 vector< double > &n_cc_proc,
                 vector< double > &n_nc_proc,
                 vector< double > &n_cc_fsi,
                 vector< double > &n_nc_fsi,
                 vector< double > &n_nue_fsi );

// -------------------------------------------------------------------------
// Make a table to compare the number of different final state particles
// and interactions in each model configuration
// -------------------------------------------------------------------------
void MakeTable( const vector<string> &model_names,
                const m_outer &n_cc_vect,
                const m_outer &n_nc_vect,
                const m_outer &n_nue_vect,
                const m_outer &n_cc_proc_vect,
                const m_outer &n_nc_proc_vect,
                const vector< string > ccinteractions,
                const vector< string > ncinteractions,
                const vector< string > nueinteractions,
                const vector< string > ccprocesses,
                const vector< string > ncprocesses,
                ostream &file );

// -------------------------------------------------------------------------
// Stacking the histograms:
// The arguments passed to the stacking function are:
//      hists            : vector of TH1Ds
//      leg_entries      : vector of strings for legend entries
//      norm             : vector of normalisation values
//      x_axis, y_axis   : axis labels
//      title            : histogram title
//      file_name        : name of the file to save the plot to
//      min, mix, n_bins : the min and max histogram x bins and the number of bins
//
// It will assign a different colour to each histogram it has to stack and 
// will set the style and min and max bin positions according to the 
// smallest min and largest max bin in the vector
// It will also take the number of bins from one of the histograms - since 
// they should all be the same, or similar
// -------------------------------------------------------------------------
void HistStacker ( vector< TH1D* >   &hists,
                   vector< string >  &leg_entries,
                   vector< double >  &norm,
                   const char* title,
                   const char* file_name,
                   const char* x_axis,
                   const char* y_axis );

// -------------------------------------------------------------------------
// Histogram stacking with statistical errors
// -------------------------------------------------------------------------
void ErrHistStacker ( vector< TH1D* >   &hists,
                   vector< string >  &leg_entries,
                   vector< double >  &norm,
                   const char* title,
                   const char* file_name,
                   const char* x_axis,
                   const char* y_axis );

// -------------------------------------------------------------------------
// Calculating the reconstructed energy of the neutrinos to compare with the 
// MC value
// -------------------------------------------------------------------------
void RecoNuE( TTree *event_tree,
              vector< double > &reco_E_CC, 
              vector< double > &reco_E_NC,
              vector< double > &MC_reco_E_CC, 
              vector< double > &MC_reco_E_NC );

// -------------------------------------------------------------------------
// Filling a section of the table
// -------------------------------------------------------------------------
void FillTableSection(const string &label,
                      const int &n_interactions,
                      const int &n_columns,
                      const vector<string> &interactions,
                      const vector<string> &mod_names,
                      const m_outer &rates,
                      ostream &ofile );


// -------------------------------------------------------------------------


#endif

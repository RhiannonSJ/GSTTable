/*
 * The functions used throughout the table production
 *
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : December 2020
 *
 */

#include "function_defs.h"

using namespace std; 

// -------------------------------------------------------------------------
//                      Normalisation function
// -------------------------------------------------------------------------
double Norm(const int nEvents, const double nPOT, const int det){
  if(det == 0)
    return (10e20/nPOT)*nEvents;
  else if(det == 1)
    return (13.2e20/nPOT)*nEvents;
  else if(det == 2)
    return (15e20/nPOT)*nEvents;
  else{
    cerr << " Error: Unknown detector enumeration, " << det 
      << "\n Options are:"
      << "\n 0 = SBND"
      << "\n 1 = MicroBooNE"
      << "\n 2 = ICARUS"
      << "\n Exiting.";
    exit(1);
  }
}

// -------------------------------------------------------------------------
//                    Make final state particles map
// -------------------------------------------------------------------------
void FSPNumbers( TTree *event_tree,
    m_map &n_fsp ){

  // Firstly, get out the trees we want to look at
  // All correspond to number of particles AFTER FSI
  //  - fspl  == 13 (muon)
  //  - fspl  == 14 (numu)
  //  - nfp   == protons
  //  - nfn   == neutrons
  //  - nfpip == pi+ } Collate these to 
  //  - nfpim == pi- } nfcpi (#final charged pions
  //  - nfpi0 == pi0
  // Set the branch addresses for these leaves
TBranch *b_fspl  = event_tree->GetBranch("fspl");
TBranch *b_nfp   = event_tree->GetBranch("nfp");
TBranch *b_nfn   = event_tree->GetBranch("nfn");
TBranch *b_nfpip = event_tree->GetBranch("nfpip");
TBranch *b_nfpim = event_tree->GetBranch("nfpim");
TBranch *b_nfpi0 = event_tree->GetBranch("nfpi0");

// Create the variables to use as counters
int nfmu   = 0;
int nfe    = 0;
int nfnumu = 0;
int nfp    = 0;
int nfn    = 0;
int nfcpi  = 0;
int nfpi0  = 0;

// Get the number of events which contain final state muons
int n_values = event_tree->GetEntries(); // Number of entries to loop over

// Loop over the leaves and calculate the reconstructed energy
for( int i = 0; i < n_values; ++i){

  // Get the current entry
  event_tree->GetEntry(i);

  // Count #final state leptons
  if( b_fspl->GetLeaf("fspl")->GetValue() == 13 ){
    ++nfmu;
  }
  else if( b_fspl->GetLeaf("fspl")->GetValue() == 14 ){
    ++nfnumu;
  }
  else if( b_fspl->GetLeaf("fspl")->GetValue() == 11 ){
    ++nfe;
  }

  // Count #final state nucleons
  if( b_nfp->GetLeaf("nfp")->GetValue() == 1 ){
    ++nfp;
  }
  if( b_nfn->GetLeaf("nfn")->GetValue() == 1 ){
    ++nfn;
  }

  // Count #final state pions
  if( b_nfpip->GetLeaf("nfpip")->GetValue() == 1 || b_nfpim->GetLeaf("nfpim")->GetValue() == 1){
    ++nfcpi;           
  }
  if( b_nfpi0->GetLeaf("nfpi0")->GetValue() == 1 ){
    ++nfpi0;           
  }
}

// Now fill the map with the values
// Do this in such a way that the key can be printed straight into a table
n_fsp.insert( pair< string, int > (" Protons ", nfp));
n_fsp.insert( pair< string, int > (" Neutrons ", nfn));
n_fsp.insert( pair< string, int > (" Muons ", nfmu));
n_fsp.insert( pair< string, int > (" Muon Neutrinos ", nfnumu));
n_fsp.insert( pair< string, int > (" Charged Pions ", nfcpi));
n_fsp.insert( pair< string, int > (" Neutral Pions ", nfpi0));
n_fsp.insert( pair< string, int > (" Electrons ", nfe));
}
// -------------------------------------------------------------------------
//                  Make final state interactions map
// -------------------------------------------------------------------------
void FSINumbers( TTree *event_tree,
    double norm,
    vector< double > &n_cc_proc,
    vector< double > &n_nc_proc,
    vector< double > &n_cc_fsi,
    vector< double > &n_nc_fsi,
    vector< double > &n_nue_fsi){
  // Firstly, get out the trees we want to look at
  // Need both cc and nc with varying number of outgoing pions
  //      - cc : charged current FSI
  //      - nc : neutral current FSI

  // Set the branch addresses for these leaves
  int n;
  int neu;
  int nuance;
  bool cc;
  bool nc;
  bool res;
  bool coh;
  bool qel;
  bool dis;
  bool mec;
  bool cha;
  int pi0;
  int pip;
  int pim;
  int k0; 
  int kp; 
  int km; 
  int p;
  int str;
  int pdgf[1000];
  double en[1000];

  event_tree->SetBranchAddress("nuance_code",    &nuance);
  event_tree->SetBranchAddress("pdgf",    &pdgf);
  event_tree->SetBranchAddress("Ef",      &en);
  event_tree->SetBranchAddress("nf",      &n);
  event_tree->SetBranchAddress("neu",     &neu);
  event_tree->SetBranchAddress("cc",      &cc);
  event_tree->SetBranchAddress("nc",      &nc);
  event_tree->SetBranchAddress("coh",     &coh);
  event_tree->SetBranchAddress("res",     &res);
  event_tree->SetBranchAddress("qel",     &qel);
  event_tree->SetBranchAddress("dis",     &dis);
  event_tree->SetBranchAddress("mec",     &mec);
  event_tree->SetBranchAddress("charm",   &cha);
  event_tree->SetBranchAddress("nfother", &str);
  event_tree->SetBranchAddress("nfpi0",   &pi0);
  event_tree->SetBranchAddress("nfpip",   &pip);
  event_tree->SetBranchAddress("nfpim",   &pim);
  event_tree->SetBranchAddress("nfk0",    &k0);
  event_tree->SetBranchAddress("nfkp",    &kp);
  event_tree->SetBranchAddress("nfkm",    &km);
  event_tree->SetBranchAddress("nfp",     &p);

  // Charged current counters
  int ncc        = 0;
  int ncc0pi     = 0;
  int ncc0pi0p   = 0;
  int ncc0pi1p   = 0;
  int ncc0pi1p50 = 0;
  int ncc0pi1p20 = 0;
  int ncc0pi2p   = 0;
  int ncc0pi2p50 = 0;
  int ncc0pi2p20 = 0;
  int ncc0pi3p   = 0;
  int ncc0pi3p50 = 0;
  int ncc0pi3p20 = 0;
  int ncc0pimp   = 0;
  int ncc0pimp50 = 0;
  int ncc0pimp20 = 0;
  int ncc1pip    = 0;
  int ncc1pim    = 0;
  int ncc1pi0    = 0;
  int ncc2pi     = 0;
  int ncc3pi     = 0;
  int nck        = 0;
  int nnk        = 0;
  int nsigpp     = 0;
  int nsigp      = 0;
  int nlamp      = 0;
  int ndmesons   = 0;
  int nmoremu    = 0;
  int nsinglek   = 0;
  int nchbaryons = 0;
  int nhyperons  = 0;
  int nlamonly   = 0;
  int nsigonly   = 0;

  int ncccoh     = 0;
  int nccres     = 0;
  int nccqel     = 0;
  int nccdis     = 0;
  int nccmec     = 0;
  int nccqecha   = 0;
  int nccqestr   = 0;
  int nccdicha   = 0;
  int nccotherpr = 0;
  int nnccoh     = 0;
  int nncres     = 0;
  int nncqel     = 0;
  int nncdis     = 0;
  int nncmec     = 0;
  int nncotherpr = 0;

  // Neutral current counters
  int nnc        = 0;
  int nnc0pi     = 0;
  int nnc1pipm   = 0;
  int nncg2pipm  = 0;
  int nncg1pi0   = 0;

  // Nue counter
  int nnue       = 0;

  // Get the number of events which contain final state muons
  int n_values = event_tree->GetEntries(); // Number of entries to loop over

  // Proton mass
  double m_proton = 0.938272; // GeV

  // Loop and count for various conditions
  for( int i = 0; i < n_values; ++i ){

    // Get the current entry
    event_tree->GetEntry(i);

    int n_p_50 = 0;
    int n_p_20 = 0;
    int nmu    = 0;

    // CC Inclusive
    if (cc){
      // Numu content only for FSI breakdown
      if(abs(neu) == 14){
        ++ncc;
        for ( int j = 0; j < n; ++ j){
          int pdg  = pdgf[j];
          double e = en[j];
          if ( pdg == 2212 ){
            double t = e - m_proton;
            if ( t > 0.05 )
              ++n_p_50;
            if ( t > 0.02 )
              ++n_p_20;
          }
          if(abs(pdg) == 13)
            ++nmu;
        }
        // Charged current
        // CC 0Pi
        if ((pi0+pip+pim == 0)){
          ++ncc0pi; 
          if (p == 0)
            ++ncc0pi0p;
          else if (p == 1)
            ++ncc0pi1p;
          else if (p == 2)
            ++ncc0pi2p;
          else if (p == 3)
            ++ncc0pi3p;
          else if (p > 3)
            ++ncc0pimp;

          // Thresholds
          if(n_p_20 != 0){
            if (n_p_20 == 1)
              ++ncc0pi1p20;
            else if (n_p_20 == 2)
              ++ncc0pi2p20;
            else if (n_p_20 == 3)
              ++ncc0pi3p20;
            else
              ++ncc0pimp20;
          }
          if(n_p_50 != 0){
            if (n_p_50 == 1)
              ++ncc0pi1p50;
            else if (n_p_50 == 2)
              ++ncc0pi2p50;
            else if (n_p_50 == 3)
              ++ncc0pi3p50;
            else
              ++ncc0pimp50;
          }
        } // CC 0pi

        // CC 1pi
        if (pip == 1 && pim == 0 && pi0 == 0)
          ++ncc1pip;
        else if (pip == 0 && pim == 1 && pi0 == 0)
          ++ncc1pim;
        else if (pip == 0 && pim == 0 && pi0 == 1)
          ++ncc1pi0;

        // CC 2pi
        else if (pip + pim + pi0 == 2)
          ++ncc2pi;

        // CC >=3pi
        else if (pip + pim + pi0 >= 3)
          ++ncc3pi;

        // CC >1 mu
        if(nmu > 1)
          ++nmoremu;

        // CC K+K-
        if (kp + km == 2)
          ++nck;

        // CC K0K0bar
        if (k0 == 2)
          ++nnk;

        // CC single kaon production
        if(kp + km + k0 == 1)
          nsinglek++;

        // For charm CCQE 
        int sigmapp_count = 0;
        int sigmap_count  = 0;
        int lambda_count  = 0;
        int lamonly_count = 0;
        int sigonly_count = 0;

        // For charm DIS
        int d_count = 0;

        for (int j = 0; j < n; ++j){
          int pdg  = pdgf[j];
          // sigma ++
          if (pdg == 4222)
            ++sigmapp_count;
          // sigma +
          else if (pdg == 4212)
            ++sigmap_count;

          // lambda +
          else if (pdg == 4122)
            ++lambda_count;

          // D
          else if (abs(pdg) == 411 || pdg == 421)
            ++d_count;
        }

        // CC Sigma++
        if (sigmapp_count == 1)
          ++nsigpp;

        // CC Sigma+
        if (sigmap_count == 1)
          ++nsigp;

        // CC Lambda+
        if (lambda_count == 1){
          ++nlamp;
        }

        if(d_count >= 1)
          ++ndmesons;

        // Hyperon production 
        int sigmapm0_count = 0;
        int lambda0_count  = 0;
        
        for (int j = 0; j < n; ++j){
          int pdg  = pdgf[j];
          // sigma +/-
          if (pdg == 3222 || pdg == 3112 || pdg == 3212){
            ++sigmapm0_count;
            if(((kp + km + k0) == 0))
              ++sigonly_count;
          }

          // lambda 0
          else if (pdg == 3122){
            ++lambda0_count;
            if(((kp + km + k0) == 0))
              ++lamonly_count;
          }
        }

        // CC hyperon
        if(sigmapm0_count > 0 || lambda0_count > 0)
          nhyperons++;
        if(lamonly_count == 1)
          nlamonly++;
        if(sigonly_count == 1)
          nsigonly++;

      } // numu
      else if(abs(neu) == 12)
        nnue++;

      // Physical processes
      if ( qel ){
        if(!cha && !str)
          ++nccqel;
        else if (cha) {
          ++nccotherpr;
          //          ++nccqecha;
        }
        else if(!cha && str > 0){
          ++nccotherpr;
          //          ++nccqestr;
        }
      }
      else if (mec)
        ++nccmec;
      else if (res)
        ++nccres;
      else if (dis){
        // if(!d_count)
        ++nccdis;
        //else if(d_count == 1)
        //  ++nccdicha;
      }
      else if (coh)
        ++ncccoh;
      else
        ++nccotherpr;
    } // cc

    // Neutral current
    else if (nc){
      // Numu content only for FSI breakdown
      if(abs(neu) == 14){
        ++nnc;
        // NC 0pi
        if (pip + pim + pi0 == 0)
          ++nnc0pi;

        // NC 1pipm
        else if ((pip + pim == 1) && pi0 == 0)
          ++nnc1pipm;

        // NC >1pi0
        else if ((pip + pim == 0) && pi0 >= 1)
          ++nncg1pi0;

        // NC >2pipm
        else if (pip + pim + pi0 >= 2)
          ++nncg2pipm;

      } // numu
      else if(abs(neu) == 12)
        nnue++;

      // Physical processes
      if ( qel ){
        if(!cha && !str)
          ++nncqel;
      }
      else if (mec)
        ++nncmec;
      else if (res)
        ++nncres;
      else if (dis)
        ++nncdis;
      else if (coh)
        ++nnccoh;
      else
        ++nncotherpr;
    } // nc
  }
  nchbaryons = nsigpp+nsigp+nlamp;

  // Now fill the map with the values
  // Do this in such a way that the key can be printed straight into a table
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi0p) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi1p) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi1p20) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi1p50) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi2p) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi2p20) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi2p50) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi3p) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi3p20) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi3p50) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pimp) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pimp20) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pimp50) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pip) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pim) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pi0) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc2pi) );
  n_cc_fsi.push_back ( TMath::Floor( norm * ncc3pi) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nmoremu) );
  //n_cc_fsi.push_back ( TMath::Floor( norm * ndmesons) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nsinglek) );
  n_cc_fsi.push_back ( TMath::Floor( norm * (nck+nnk)) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nhyperons) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nlamonly) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nsigonly) );
  n_cc_fsi.push_back ( TMath::Floor( norm * nchbaryons) );

  n_cc_proc.push_back ( TMath::Floor( norm * nccqel) );
  n_cc_proc.push_back ( TMath::Floor( norm * nccmec) );
  n_cc_proc.push_back ( TMath::Floor( norm * nccres) );
  n_cc_proc.push_back ( TMath::Floor( norm * nccdis) );
  n_cc_proc.push_back ( TMath::Floor( norm * ncccoh) );
//  n_cc_proc.push_back ( TMath::Floor( norm * nccqecha) );
//  n_cc_proc.push_back ( TMath::Floor( norm * nccqestr) );
  n_cc_proc.push_back ( TMath::Floor( norm * nccotherpr) );
  //n_cc_proc.push_back ( TMath::Floor( norm * nccdicha) );
  n_nc_proc.push_back ( TMath::Floor( norm * nncqel) );
  n_nc_proc.push_back ( TMath::Floor( norm * nncmec) );
  n_nc_proc.push_back ( TMath::Floor( norm * nncres) );
  n_nc_proc.push_back ( TMath::Floor( norm * nncdis) );
  n_nc_proc.push_back ( TMath::Floor( norm * nnccoh) );
  n_nc_proc.push_back ( TMath::Floor( norm * nncotherpr) );

  n_nc_fsi.push_back ( TMath::Floor( norm * nnc) );
  n_nc_fsi.push_back ( TMath::Floor( norm * nnc0pi) );
  n_nc_fsi.push_back ( TMath::Floor( norm * nnc1pipm) );
  n_nc_fsi.push_back ( TMath::Floor( norm * nncg2pipm) );
  n_nc_fsi.push_back ( TMath::Floor( norm * nncg1pi0) );
  //
  n_nue_fsi.push_back( TMath::Floor( norm * nnue) );

}

// -------------------------------------------------------------------------
//                      Make final state tables
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
                ostream &file ){

  // Iterators 
  typedef m_outer::const_iterator map_it;

  // Make a table from an input map - of - maps to print nice things
  // Number of columns and rows to be made
  int n_models, n_ccinteractions, n_ncinteractions, n_nueinteractions, n_ccprocesses, n_ncprocesses, n_cols;

  n_models          = model_names.size();
  n_ccinteractions  = ccinteractions.size();
  n_ncinteractions  = ncinteractions.size();
  n_nueinteractions = nueinteractions.size();
  n_ccprocesses     = ccprocesses.size();
  n_ncprocesses     = ncprocesses.size();

  // If we are only looking at a single model, add a column with the statistical uncertainty
  if(n_models == 1) n_cols = 2;
  else n_cols = n_models;

                    
  // Begin the tabular environment in the LaTeX output file for n_predictions columns
  file << "\\begin{longtable}{ m{4cm} * {" << n_cols << "}{ >{\\centering\\arraybackslash}m{2.8cm} } }" << endl;
  file << "\\hline" << endl;

  // Fill the first line with "Configurations" for all the right-hand columns
  file << " \\multirow{2}{*}{\\textbf{ Hadronic Final State }} & \\multicolumn{ " << n_cols << " }{ c }{ \\textbf{ Model Configurations } } \\\\" << endl;


  // Fill a line of the table with the prediction names
  file << " & " ;
  for( int i = 0; i < n_cols-1; ++i ){
    file << " \\textbf{ " << model_names.at(i) << " } & ";
  }
  if(n_cols != n_models) // Then we want the statistical uncertainty as the final column
    file << " \\textbf{ Stat. Err. } \\\\" << endl;
  else
    file << " \\textbf{ " << model_names.at(n_models-1) << " } \\\\" << endl;

  file << " \\hline " << endl;

  // -------------------------------------------------------------------
  //      Loop over the outer maps and fill the rest of the table                                   
  // -------------------------------------------------------------------
  FillTableSection("$\\nu_{\\mu}$ \\textit{ Charged Current }", n_ccinteractions, n_cols, ccinteractions, model_names, n_cc_vect, file);
  file << " \\hdashline " << endl;
  FillTableSection("$\\nu_{\\mu}$ \\textit{ Neutral Current }", n_ncinteractions, n_cols, ncinteractions, model_names, n_nc_vect, file);
  file << " \\hdashline " << endl;
  FillTableSection("$\\nu_{e}$", n_nueinteractions, n_cols, nueinteractions, model_names, n_nue_vect, file);

  file << " \\hline " << endl;
  file << " \\textbf{Physical Process } &" << endl;
  for ( int i = 0; i < n_cols - 1; ++i ){
    file << " & " << endl;
  }
  file << " \\\\ " << endl;
  file << " \\hline " << endl;

  FillTableSection("\\textit{ Charged Current }", n_ccprocesses, n_cols, ccprocesses, model_names, n_cc_proc_vect, file);
  file << " \\hdashline " << endl;
  FillTableSection("\\textit{ Neutral Current }", n_ncprocesses, n_cols, ncprocesses, model_names, n_nc_proc_vect, file);

  file << " \\hline " << endl;

  // -------------------------------------------------------------------
  //                    End the table                                   
  // -------------------------------------------------------------------

  // End the tabular environment
  file << "\\end{longtable}" << endl;


}

// -------------------------------------------------------------------------
//                      hist stacking function
// -------------------------------------------------------------------------
void HistStacker ( vector< TH1D* >   &hists,
    vector< string >  &leg_entries,
    vector< double >  &norm, 
    const char* title,
    const char* file_name,
    const char* x_axis,
    const char* y_axis ){

  // The Canvas, legend and empty histogram to print the title and axes labels
  // Canvas
  TCanvas *c   = new TCanvas( "c", title, 800, 600 );
  TLegend *leg = new TLegend( 0.68, 0.68, 0.88, 0.88 );

  // Get the sizes of the vectors
  int n_hists, n_leg, n_norm;
  n_hists = hists.size();
  n_leg   = leg_entries.size();
  n_norm  = norm.size();

  // Check they are the same, or return an error
  if ( n_hists != n_leg || n_hists != n_norm || n_leg != n_norm ) {
    cerr << " The vectors should all have the same number of entries " << endl;
    exit(1);
  }

  // Loop over the histograms 
  for ( int i = 0; i < n_hists; ++i ) {

    // Set their line colours, styles and widths
    if ( i < 4 || i == 5 ) {
      hists[i]->SetLineColor( i + 1 );
      hists[i]->SetLineStyle( 1 );
      hists[i]->SetLineWidth( 1 );
      hists[i]->SetMarkerColor( i + 1 );
    }
    else if ( i == 4 || i == 6 || i == 7 ){
      hists[i]->SetLineColor( i + 3 );
      hists[i]->SetLineStyle( 1 );
      hists[i]->SetLineWidth( 1 );
      hists[i]->SetMarkerColor( i + 3 );       
    }

    // Normalise the histogram with the associated normalisation
    hists[i]->Scale( norm[i] );
  }

  for ( int i = 0; i < n_hists; ++i ) {

    if( i != 4){
      // Add in the legend entry
      leg->AddEntry( hists[i], leg_entries[i].c_str(), "l" );
    }
  }

  // Variables for calculating the maximum y value
  int i_max;
  double max = -1000; 
  double max_y;

  // Loop over the histograms and find which has the highest y-val
  // Print this first
  for ( int i = 0; i < n_hists; ++i ) {
    if ( hists[i]->GetMaximum() > max){
      max = hists[i]->GetMaximum();
    }
  }
  max_y = max + 0.1*max;

  // Draw the histograms
  hists[0]->Draw();
  for ( int i = 1; i < n_hists; ++i ) {        

    // For now, don't draw G17_01b
    if( i != 4 ){
      // Draw the histograms 
      hists[i]->Draw( "same" );
    }

  }

  hists[0]->GetXaxis()->SetTitle(x_axis);
  hists[0]->GetYaxis()->SetTitle(y_axis);
  hists[0]->SetAxisRange(0,max_y, "Y");
  hists[0]->SetTitleOffset(1.5, "Y");    
  hists[0]->SetStats(kFALSE);

  leg->Draw();
  c->SaveAs(file_name);

  delete c;
  delete leg;

}

// -------------------------------------------------------------------------
//             hist stacking function with statistical errors
// -------------------------------------------------------------------------
void ErrHistStacker ( vector< TH1D* >   &hists,
    vector< string >  &leg_entries,
    vector< double >  &norm, 
    const char* title,
    const char* file_name,
    const char* x_axis,
    const char* y_axis ){

  // The Canvas, legend and empty histogram to print the title and axes labels
  // Canvas
  TCanvas *c   = new TCanvas( "c", title, 600, 600 );
  TLegend *leg = new TLegend( 0.68, 0.68, 0.88, 0.88 );

  // Get the sizes of the vectors
  int n_hists, n_leg, n_norm;
  n_hists = hists.size();
  n_leg   = leg_entries.size();
  n_norm  = norm.size();

  // Check they are the same, or return an error
  if ( n_hists != n_leg || n_hists != n_norm || n_leg != n_norm ) {
    cerr << " The vectors should all have the same number of entries " << endl;
    exit(1);
  }

  // Loop over the histograms 
  for ( int i = 0; i < n_hists; ++i ) {

    // Set their line colours, styles and widths
    if ( i < 4 || i == 5 ) {
      hists[i]->SetLineColor( i + 1 );
      hists[i]->SetLineStyle( 1 );
      hists[i]->SetLineWidth( 1 );
      hists[i]->SetMarkerColor( i + 1 );
    }
    else if ( i == 4 || i == 6 || i == 7 ){
      hists[i]->SetLineColor( i + 3 );
      hists[i]->SetLineStyle( 1 );
      hists[i]->SetLineWidth( 1 );
      hists[i]->SetMarkerColor( i + 3 );       
    }

    // Normalise the histogram with the associated normalisation
    hists[i]->Sumw2();
    hists[i]->Scale( norm[i] );
  }

  for ( int i = 0; i < n_hists; ++i ) {

    if( i != 4){
      // Add in the legend entry
      leg->AddEntry( hists[i], leg_entries[i].c_str(), "l" );
    }
  }

  // Variables for calculating the maximum y value
  int i_max;
  double max = -1000; 
  double max_y;

  // Loop over the histograms and find which has the highest y-val
  // Print this first
  for ( int i = 0; i < n_hists; ++i ) {
    if ( hists[i]->GetMaximum() > max){
      max = hists[i]->GetMaximum();
    }
  }
  max_y = max + 0.1*max;

  /*
  // ----------------------------------------------------------------------
  //                          Statistical errors
  // ----------------------------------------------------------------------
  // Loop over histograms
  for ( int i = 0; i < n_hists; ++i ){
  // Loop over bins, get value, calculate sqrt(N)
  // Vector for each hist
  int n_bins = 0;
  n_bins = hists[i]->GetNbinsX();

  for ( int j = 0; j < n_bins; ++j ){
  double n_events = 0;
  n_events = hists[i]->GetBinContent(j);   

  double err_val = 0;
  err_val = sqrt(n_events);
  hists[i]->SetBinError(j,err_val);
  }   

  }
  */
  // ----------------------------------------------------------------------
  //                                Draw
  // ----------------------------------------------------------------------


  gStyle->SetEndErrorSize(3);
  gStyle->SetErrorX(0);

  // Draw the histograms
  //hists[0]->Draw("e1");
  hists[0]->Draw("e1");
  hists[0]->Draw("same chist");
  for ( int i = 1; i < n_hists; ++i ) {        

    // For now, don't draw G17_01b
    if( i != 4 ){
      // Draw the histograms 
      hists[i]->Draw( "e1same" );
      hists[i]->Draw( "same chist" );
    }

  }

  hists[0]->GetXaxis()->SetTitle(x_axis);
  hists[0]->GetYaxis()->SetTitle(y_axis);
  hists[0]->SetAxisRange(0,max_y, "Y");
  hists[0]->SetTitleOffset(1.5, "Y");    
  hists[0]->SetStats(kFALSE);

  leg->Draw();
  c->SaveAs(file_name);

  delete c;
  delete leg;

}

// -------------------------------------------------------------------------
//                    reconstructed energy calculation
// -------------------------------------------------------------------------
void RecoNuE( TTree *event_tree,
    vector< double > &reco_E_CC, 
    vector< double > &reco_E_NC,
    vector< double > &MC_reco_E_CC, 
    vector< double > &MC_reco_E_NC ){

  // Get the branches to calculate reconstructed energy and MC energy
  TBranch *b_mu_e  = event_tree->GetBranch("El");
  TBranch *b_nu_e  = event_tree->GetBranch("Ev");
  TBranch *b_mu_p  = event_tree->GetBranch("pl");
  TBranch *b_theta = event_tree->GetBranch("cthl");
  TBranch *b_nfpi0 = event_tree->GetBranch("nfpi0");
  TBranch *b_nfpip = event_tree->GetBranch("nfpip");
  TBranch *b_nfpim = event_tree->GetBranch("nfpim");
  TBranch *b_cc    = event_tree->GetBranch("cc");
  TBranch *b_nc    = event_tree->GetBranch("nc");

  // The variables from the branches and get the leaves
  double m_n   = 0.93828;   // Nucleon mass, GeV
  double m_mu  = 0.10566;   // Muon mass, GeV

  int n_values = event_tree->GetEntries(); // Number of entries to loop over

  // Loop over the leaves and calculate the reconstructed energy
  for( int i = 0; i < n_values; ++i){

    event_tree->GetEntry(i);

    double reco, reco_mc, e, p, cth;

    // For CC0pi
    if( b_cc->GetLeaf("cc")->GetValue() != 0 
        && b_nfpip->GetLeaf("nfpip")->GetValue()
        + b_nfpim->GetLeaf("nfpim")->GetValue()
        + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0 ){

      // Get the values needed
      e   = b_mu_e->GetLeaf("El")->GetValue();
      p   = b_mu_p->GetLeaf("pl")->GetValue();
      cth = b_theta->GetLeaf("cthl")->GetValue(); 

      reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( e - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 

      reco_mc = reco - double(b_nu_e->GetLeaf("Ev")->GetValue());

      // Make the vectors of reconstructed and reconstructed-MC energy for CC0pi
      MC_reco_E_CC.push_back(reco_mc);
      reco_E_CC.push_back(reco);

    }
    // For NC0pi
    else if( b_nc->GetLeaf("nc")->GetValue() != 0 
        && b_nfpip->GetLeaf("nfpip")->GetValue()
        + b_nfpim->GetLeaf("nfpim")->GetValue()
        + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0 ){

      // Get the values
      e   = b_mu_e->GetLeaf("El")->GetValue();
      p   = b_mu_p->GetLeaf("pl")->GetValue();
      cth = b_theta->GetLeaf("cthl")->GetValue(); 

      reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( e - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 

      reco_mc = reco - double(b_nu_e->GetLeaf("Ev")->GetValue());

      // Make the vectors of reconstructed and reconstructed-MC energy for CC0pi
      MC_reco_E_NC.push_back(reco_mc);
      reco_E_NC.push_back(reco);

    }
  }
}
// -------------------------------------------------------------------------
//                    Filling a table section
// -------------------------------------------------------------------------
void FillTableSection(const string &label,
                      const int &n_interactions,
                      const int &n_columns,
                      const vector<string> &interactions,
                      const vector<string> &mod_names,
                      const m_outer &rates,
                      ostream &ofile){

  ofile << " \\multicolumn{ " << n_columns + 1 << " }{  c  }{ " << label << " } \\\\ " << endl;
  ofile << " \\hdashline " << endl;

  const int n_mods = mod_names.size();
  for( int i = 0; i < n_interactions; ++i ){

    ofile << interactions.at(i) << " & ";

    for(int n = 0; n < n_columns - 1; ++n){
      std::string m = mod_names.at(n);

      ofile << std::fixed << setprecision(0) <<  "\\num{ " << rates.at(m).at(i) << "} & ";
    }
    if(n_columns != n_mods){
      std::string m = mod_names.at(0);
      float stat_err = 100./static_cast<float>(std::sqrt(rates.at(m).at(i)));
      ofile << std::fixed << setprecision(2) << stat_err << "\\% \\\\ " << endl;
    }
    else{
      std::string m = mod_names.at(n_mods-1);
      ofile << std::fixed << setprecision(0) << "\\num{ " <<  rates.at(m).at(i) << "} \\\\ " << endl;
    }
  } 

}
void function_defs(){}

/*
 *
 * A macro to draw a LATEX table containing event rate multiplicities 
 * from vaious GENIE models
 *
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : December 2020
 *
*/

#include "function_defs.h"

void table(const std::string &fileList, const int &det=0, const std::string &outputFile="outputTable.tex", const bool avg=false) {
    
  // -------------------------------------------------------------------------
  //        Loop over fileList and calculate everything to fill table
  // -------------------------------------------------------------------------

  //Construct the list of model and file names from the input
  ifstream files(fileList.c_str());
  std::string line;
  std::vector<std::string> lines, modelNames, fileNames;
  while(getline(files,line)){
    lines.push_back(line);
  }

  // Now extract the model and file names from the list
  for(const std::string &l : lines){
    int pos0=0, pos1;
    std::string modelName="", fileName="";
    while(l[pos0] != ' '){
      modelName += l[pos0];
      pos0++;
    }
    for(pos1 = pos0; l[pos1] == ' '; ++pos1);
    while(l[pos1] != ' '){
      fileName += l[pos1];
      pos1++;
    }
    modelNames.push_back(modelName);
    fileNames.push_back(fileName);
  }
  assert(modelNames.size() == fileNames.size());

  // Initial and final CC and NC maps
  m_outer cc_proc_model_ints;
  m_outer nc_proc_model_ints;
  m_outer cc_reco_model_ints;
  m_outer nc_reco_model_ints;
  m_outer nue_reco_model_ints;

  // Loop over the files and get the average POT
  double avg_pot = 0.;
  for(unsigned int i = 0; i < modelNames.size(); ++i){
    TFile f(fileNames.at(i).c_str(),"READ");
    if (f.IsZombie()) {
      cout << " Error opening file " << fileNames.at(i) << endl;
      exit(-1);
    }
    // Get the trees we want from the root files
    TTree *gst   = static_cast<TTree*>(f.Get("gst"));
    double pot   = gst->GetWeight();
    avg_pot     += pot;
  }
  if(avg) avg_pot /= static_cast<double>(modelNames.size());

  // Now fill the relevant quantities for the table
  for(unsigned int i = 0; i < modelNames.size(); ++i){
    TFile f(fileNames.at(i).c_str(),"READ");
    if (f.IsZombie()) {
      cout << " Error opening file " << fileNames.at(i) << endl;
      exit(-1);
    }
    else{
      cout << " " << modelNames.at(i) << " file is open " << endl;
    }
    // Get the trees we want from the root files
    TTree *gst   = static_cast<TTree*>(f.Get("gst"));
    int    nev = gst->GetEntries();

    std::cout << " Average POT across all TTrees: " << avg_pot << " corresponding to " << nev << " events" << std::endl; 

    // -------------------------------------------------------------------------
    //                          Get the normalisations
    // -------------------------------------------------------------------------

    double norm     = Norm(1,avg_pot,det);
    double normRate = Norm(nev,avg_pot,det);
    std::cout << " Total events in the active volume from " << modelNames.at(i) << " model is: " << normRate << std::endl;
    std::cout << " Normalisation to apply to each event from " << modelNames.at(i) << " model is: " << norm << std::endl;

    vector< double > n_cc_fsi, n_nc_fsi, n_nue_fsi, n_cc_proc, n_nc_proc;
    FSINumbers(gst, norm, n_cc_proc, n_nc_proc, n_cc_fsi, n_nc_fsi, n_nue_fsi);
  
    cc_reco_model_ints .emplace(modelNames.at(i), n_cc_fsi  ); 
    nc_reco_model_ints .emplace(modelNames.at(i), n_nc_fsi  ); 
    nue_reco_model_ints.emplace(modelNames.at(i), n_nue_fsi );
    cc_proc_model_ints .emplace(modelNames.at(i), n_cc_proc );
    nc_proc_model_ints .emplace(modelNames.at(i), n_nc_proc );

  }
  std::cout << "-----------------------------------------------------" << std::endl;
  // Make vector of final hadronic states
  // Using LaTeX syntax
  vector< string > CCFHS, NCFHS, NUEFHS, CCPP, NCPP;

  CCFHS.push_back( "Inclusive" );
  CCFHS.push_back( "\\hspace{.3cm} 0$ \\pi$ + X" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 0p" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 1p" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 1p ($ > $ 20 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 1p ($ > $ 50 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 2p" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 2p ($ > $ 20 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 2p ($ > $ 50 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 3p" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 3p ($ > $ 20 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi $ 3p ($ > $ 50 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi > $ 3p" );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi > $ 3p ($ > $ 20 MeV) " );
  CCFHS.push_back( "\\hspace{.6cm} 0$ \\pi > $ 3p ($ > $ 50 MeV) " );    
  CCFHS.push_back( "\\hspace{.3cm} 1$ \\pi^+ $ + X" );
  CCFHS.push_back( "\\hspace{.3cm} 1$ \\pi^- $ + X" );
  CCFHS.push_back( "\\hspace{.3cm} 1$ \\pi^0 $ + X" );
  CCFHS.push_back( "\\hspace{.3cm} 2$ \\pi $ + X" );
  CCFHS.push_back( "\\hspace{.3cm} $ \\geqslant 3 \\pi $ + X" );
  CCFHS.push_back( "\\hspace{.3cm} $ > 1 \\mu $ + X" );
  //CCFHS.push_back( "\\hspace{.3cm} $ D^{\\pm}, D^{0} $ (+ X)" );
  CCFHS.push_back( "\\hspace{.3cm} $ K^{\\pm}, K^{0} $ (+ X)" );
  CCFHS.push_back( "\\hspace{.3cm} $ K^{+} K^{-}, K^{0} \\bar{K}^{0} $ (+ X)" );
  CCFHS.push_back( "\\hspace{.3cm} $ \\Sigma^{\\pm}, \\Sigma^{0}, \\Lambda^{0} $ (+ X)" );
  CCFHS.push_back( "\\hspace{.3cm} $ \\Lambda^{0} $" );
  CCFHS.push_back( "\\hspace{.3cm} $ \\Sigma^{0}, \\Sigma^{\\pm} $" );
  CCFHS.push_back( "\\hspace{.3cm} $ \\Sigma^{++}_{c}, \\Sigma^{+}_{c}, \\Lambda^{+}_{c} $ (+ X)" );

  NCFHS.push_back( "Inclusive" );
  NCFHS.push_back( "\\hspace{.3cm} 0$ \\pi $" );
  NCFHS.push_back( "\\hspace{.3cm} 1$ \\pi^{\\pm} $ + X" );
  NCFHS.push_back( "\\hspace{.3cm}  $ \\geqslant 2 \\pi^{\\pm} $ + X" );
  NCFHS.push_back( "\\hspace{.3cm}  $ \\geqslant 1 \\pi^0 $ + X" );
  
  NUEFHS.push_back( "Inclusive" );

  CCPP.push_back( "QE" );
  CCPP.push_back( "MEC" );
  CCPP.push_back( "RES" );
  CCPP.push_back( "DIS" );
  CCPP.push_back( "Coherent" );
  CCPP.push_back( "Other" );

  NCPP.push_back( "QE" );
  NCPP.push_back( "MEC" );
  NCPP.push_back( "RES" );
  NCPP.push_back( "DIS" );
  NCPP.push_back( "Coherent" );
  NCPP.push_back( "Other" );

  ofstream file_reco;
  file_reco.open( outputFile.c_str() );
  MakeTable(modelNames, cc_reco_model_ints, nc_reco_model_ints, nue_reco_model_ints, cc_proc_model_ints, nc_proc_model_ints, CCFHS, NCFHS, NUEFHS, CCPP, NCPP, file_reco );

}

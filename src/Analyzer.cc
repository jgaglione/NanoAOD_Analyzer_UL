#include "Analyzer.h"
#include "Compression.h"
#include <regex>
#include <sstream>
#include <cmath>
#include <map>
#include <typeinfo>
//// Used to convert Enums to integers
#define ival(x) static_cast<int>(x)
//// BIG_NUM = sqrt(sizeof(int)) so can use diparticle convention of
//// index = BIG_NUM * i1 + i2
//// This insures easy way to extract indices
//// Needs to be changed if go to size_t instead (if want to play safe
#define BIG_NUM 46340

///// Macros defined to shorten code.  Made since lines used A LOT and repeative.  May change to inlines
///// if tests show no loss in speed
#define histAddVal2(val1, val2, name) ihisto.addVal(val1, val2, group, max, name, wgt);
#define histAddVal(val, name) ihisto.addVal(val, group, max, name, wgt);
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef std::vector<int>::iterator vec_iter;



//////////////////////////////////////////////////////////////////
///////////////////CONSTANTS DEFINITONS///////////////////////////
//////////////////////////////////////////////////////////////////

//Filespace that has all of the .in files
const std::string PUSPACE = "Pileup/";


//////////PUBLIC FUNCTIONS////////////////////

const std::vector<CUTS> Analyzer::genCuts = {
  CUTS::eGTau, CUTS::eGNuTau, CUTS::eGTop,
  CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ,
  CUTS::eGW, CUTS::eGHiggs, CUTS::eGJet, CUTS::eGBJet, CUTS::eGHadTau, CUTS::eGMatchedHadTau
};

const std::vector<CUTS> Analyzer::jetCuts = {
  CUTS::eRJet1,  CUTS::eRJet2,   CUTS::eRCenJet,
  CUTS::eR1stJet, CUTS::eR2ndJet, CUTS::eRBJet
};

const std::vector<CUTS> Analyzer::nonParticleCuts = {
  CUTS::eRVertex,CUTS::eRTrig1, CUTS::eRTrig2,
};
//01.16.19
const std::unordered_map<std::string, CUTS> Analyzer::cut_num = {
  {"NGenTau", CUTS::eGTau},                             {"NGenTop", CUTS::eGTop},
  {"NGenElectron", CUTS::eGElec},                       {"NGenMuon", CUTS::eGMuon},
  {"NGenZ", CUTS::eGZ},                                 {"NGenW", CUTS::eGW},
  {"NGenHiggs", CUTS::eGHiggs},                         {"NGenJet", CUTS::eGJet},
  {"NGenBJet", CUTS::eGBJet},                           {"NGenHadTau", CUTS::eGHadTau},
  {"NMatchedGenHadTau", CUTS::eGMatchedHadTau},
  {"NRecoMuon1", CUTS::eRMuon1},                        {"NRecoMuon2", CUTS::eRMuon2},
  {"NRecoElectron1", CUTS::eRElec1},                    {"NRecoElectron2",CUTS::eRElec2},
  {"NRecoTau1", CUTS::eRTau1},                          {"NRecoTau2", CUTS::eRTau2},
  {"NRecoJet1", CUTS::eRJet1},                          {"NRecoJet2", CUTS::eRJet2},
  {"NRecoCentralJet", CUTS::eRCenJet},                  {"NRecoBJet", CUTS::eRBJet},
  {"NRecoTriggers1", CUTS::eRTrig1},                    {"NRecoTriggers2", CUTS::eRTrig2},
  {"NRecoFirstLeadingJet", CUTS::eR1stJet},             {"NRecoSecondLeadingJet", CUTS::eR2ndJet},
  {"NDiMuonCombinations", CUTS::eDiMuon},               {"NDiElectronCombinations", CUTS::eDiElec},
  {"NDiTauCombinations", CUTS::eDiTau},                 {"NDiJetCombinations", CUTS::eDiJet},
  {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1},         {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2},
  {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1},         {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
  {"NElectron1Tau1Combinations", CUTS::eElec1Tau1},     {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
  {"NElectron2Tau1Combinations", CUTS::eElec2Tau1},     {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
  {"NMuon1Electron1Combinations", CUTS::eMuon1Elec1},   {"NMuon1Electron2Combinations", CUTS::eMuon1Elec2},
  {"NMuon2Electron1Combinations", CUTS::eMuon2Elec1},   {"NMuon2Electron2Combinations", CUTS::eMuon2Elec2},
  {"NElectron1Jet1Combinations", CUTS::eElec1Jet1},     {"NElectron1Jet2Combinations", CUTS::eElec1Jet2},
  {"NElectron2Jet1Combinations", CUTS::eElec2Jet1},     {"NElectron2Jet2Combinations", CUTS::eElec2Jet2},
  {"NLeadJetCombinations", CUTS::eSusyCom},             {"METCut", CUTS::eMET},
  {"NRecoWJet", CUTS::eRWjet},                          {"NRecoVertex", CUTS::eRVertex},
  {"EENoiseVeto", CUTS::eREENoiseVeto}
};

const std::map<PType, float> leptonmasses = {
     {PType::Electron, 0.000511},
     {PType::Muon, 0.1056},
     {PType::Tau, 1.777}
};

		// SUSY SOFT IDs: Apply a veto to this jet to which we further apply pt and eta cuts and b-tagging cuts (Loose WPs from DeepCSV in all years, custom in 2016):
const std::map<std::string, float> susysoftlepidloosebtagdiscr = {
	{"2016", 0.4}, {"2017", 0.1522}, {"2018", 0.1241}
};
//////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS////////////////////
//////////////////////////////////////////////////////

///Constructor
Analyzer::Analyzer(std::vector<std::string> infiles, std::string outfile, bool setCR, std::string configFolder, std::string year) : goodParts(getArray()), genName_regex(".*([A-Z][^[:space:]]+)"){
  std::cout << "setup start" << std::endl;

  routfile = TFile::Open(outfile.c_str(), "RECREATE", outfile.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));

  if(distats["Run"].bfind("AddMetaData")) add_metadata(infiles);

  BOOM= new TChain("Events");
  infoFile=0;

  for(std::string infile: infiles){
    BOOM->AddFile(infile.c_str());
  }


  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  srand(0);

  filespace=configFolder;//"PartDet";
  filespace+="/";

  setupGeneral(year);

  CalculatePUSystematics = distats["Run"].bfind("CalculatePUSystematics");
  // New variable to do special PU weight calculation (2017)
  specialPUcalculation = distats["Run"].bfind("SpecialMCPUCalculation");

  // If data, always set specialcalculation to false
  if(isData){
    specialPUcalculation = false;
    std::cout << "This is Data!! Setting SpecialMCPUCalculation to False." << std::endl;
  }

  // New! Initialize pileup information and take care of exceptions if histograms not found.
  initializePileupInfo(specialPUcalculation, outfile);

  if(distats["Run"].bfind("ApplyNPVWeight")) initializeNPVWeights(year);

  syst_names.push_back("orig");
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> tmp;
  syst_parts.push_back(tmp);
  if(!isData && distats["Systematics"].bfind("useSystematics")) {
    for(auto systname : distats["Systematics"].bset) {
      if( systname == "useSystematics")
        doSystematics= true;
      else {
        syst_names.push_back(systname);
        syst_parts.push_back(getArray());
      }
    }
  }else {
    doSystematics=false;
  }

  _Electron = new Electron(BOOM, filespace + "Electron_info.in", syst_names, year);

  _Muon     = new Muon(BOOM, filespace + "Muon_info.in", syst_names, year);

  _Tau      = new Taus(BOOM, filespace + "Tau_info.in", syst_names, year);

  _Jet      = new Jet(BOOM, filespace + "Jet_info.in", syst_names, year);

  _FatJet   = new FatJet(BOOM, filespace + "FatJet_info.in", syst_names, year);

  _Photon   = new Photon(BOOM, filespace + "Photon_info.in", syst_names, year);

  if(year.compare("2017") == 0){
    _MET      = new Met(BOOM, "METFixEE2017" , syst_names, distats["Run"].dmap.at("MT2Mass"));
  }
  else{
    _MET      = new Met(BOOM, "MET" , syst_names, distats["Run"].dmap.at("MT2Mass"));
  }


  // B-tagging scale factor stuff
  setupBJetSFInfo(_Jet->pstats["BJet"], year);

  // Tau scale factors stuff
  setupTauIDSFsInfo(_Tau->pstats["TauID"].smap.at("TauIDAlgorithm"), year, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("TauSFforEmbeddedSamples"));
  setupTauResSFsInfo(distats["Run"].bfind("ApplyETauFakeRateESSF"));

  if(!isData) {
    std::cout<<"This is MC if not, change the flag!"<<std::endl;
    _Gen = new Generated(BOOM, filespace + "Gen_info.in", syst_names);
    _GenHadTau = new GenHadronicTaus(BOOM, filespace + "Gen_info.in", syst_names);
    _GenJet = new GenJets(BOOM, filespace + "Gen_info.in", syst_names);
     allParticles= {_Gen,_GenHadTau,_GenJet,_Electron,_Muon,_Tau,_Jet,_FatJet,_Photon};
  } else {
    std::cout<<"This is Data if not, change the flag!"<<std::endl;
    allParticles= {_Electron,_Muon,_Tau,_Jet,_FatJet};
  }

  particleCutMap[CUTS::eGElec]=_Electron;
  particleCutMap[CUTS::eGMuon]=_Muon;
  particleCutMap[CUTS::eGTau]=_Tau;

  std::vector<std::string> cr_variables;
  if(setCR) {
    char buf[64];
    read_info(filespace + "Control_Regions.in");
    crbins = pow(2.0, distats["Control_Region"].dmap.size());
    for(auto maper: distats["Control_Region"].dmap) {
      cr_variables.push_back(maper.first);
      sprintf(buf, "%.*G", 16, maper.second);
      cr_variables.push_back(buf);
    }
    if(isData) {
      if(distats["Control_Region"].smap.find("SR") == distats["Control_Region"].smap.end()) {
        std::cout << "Using Control Regions with data, but no signal region specified can lead to accidentially unblinding a study  before it should be.  Please specify a SR in the file PartDet/Control_Region.in" << std::endl;
        exit(1);
      } else if(distats["Control_Region"].smap.at("SR").length() != distats["Control_Region"].dmap.size()) {
        std::cout << "Signal Region specified incorrectly: check signal region variable to make sure the number of variables matches the number of signs in SR" << std::endl;
        exit(1);
      }
      int factor = 1;
      SignalRegion = 0;
      for(auto gtltSign: distats["Control_Region"].smap["SR"]) {
        if(gtltSign == '>') SignalRegion += factor;
        factor *= 2;
      }
      if(distats["Control_Region"].smap.find("Unblind") != distats["Control_Region"].smap.end()) {

        blinded = distats["Control_Region"].smap["Unblind"] == "false";
        std::cout << "we have " << blinded << std::endl;
      }
    }
  }
  //we update the root file if it exist so now we have to delete it:
  //std::remove(outfile.c_str()); // delete file
  histo = Histogramer(1, filespace+"Hist_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables);
  if(doSystematics)
    syst_histo=Histogramer(1, filespace+"Hist_syst_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables,syst_names);

  systematics = Systematics(distats);

  setupJetCorrections(year, outfile);

  if(setCR) {
    cuts_per.resize(histo.get_folders()->size());
    cuts_cumul.resize(histo.get_folders()->size());
  } else {
    cuts_per.resize(histo.get_cuts()->size());
    cuts_cumul.resize(histo.get_cuts()->size());
  }

  create_fillInfo();
  for(auto maper: distats["Control_Region"].dmap) {

    setupCR(maper.first, maper.second);
  }
  // check if we need to make gen level cuts to cross clean the samples:
  for(auto iselect : gen_selection){
    if(iselect.second){
      std::cout<<"Waning: The selection "<< iselect.first<< " is active!"<<std::endl;
    }
  }

  if(distats["Run"].bfind("DiscrByGenDileptonMass")){
    // If the gen-dilepton mass filter is on, check that the sample is DY. Otherwise, the filter won't be applied.
    isZsample = infiles[0].find("DY") != std::string::npos;
    // if(isZsample){std::cout << "This is a DY sample " << std::endl;}
    // else{std::cout << "Not a DY sample " << std::endl;}
  }

  // Check if this is a W/Z+jets sample for the purposes of applying Z-pt corrections.
  if((infiles[0].find("WJets") != std::string::npos) || (infiles[0].find("DY") != std::string::npos) || (infiles[0].find("EWKWPlus") != std::string::npos)
    || (infiles[0].find("EWKWMinus") != std::string::npos) || (infiles[0].find("Z2Jets") != std::string::npos)){
    isVSample = true;
  }
  else{
    isVSample = false;
  }

  initializeWkfactor(infiles);
  setCutNeeds();

  std::cout << "setup complete" << std::endl << std::endl;
  start = std::chrono::system_clock::now();
}

void Analyzer::add_metadata(std::vector<std::string> infiles){
  std::cout << "------------------------------------------------------------ " << std::endl;
  std::cout << "Copying minimal original trees from input files:"<<std::endl;

  // Define all the variables needed for this function
  TFile* rfile;
  TTree* keytree;

  // Loop over the list of input files.
  for( std::string infile: infiles){
    std::cout << "File: " << infile << std::endl;

    // Open the input file
    rfile = TFile::Open(infile.c_str());
    routfile->cd();
    // Loop over all key stored in the current input file
    std::cout << "Processing key: " << std::endl;
    for(const auto&& inkey : *rfile->GetListOfKeys()){
      std::string keyname = inkey->GetName();
      std::cout << "\t" << keyname << std::endl;

      if(keyname == "Events"){
        keytree= ((TTree*) rfile->Get(keyname.c_str()));     // Get the tree from file
        keytree->SetBranchStatus("*",0);                     // Disable all branches
        keytree->SetBranchStatus("run",1);                   // Enable only the branch named run
        originalTrees[keyname] = keytree->CopyTree("1","",1); // Add this tree to the original trees map. No selection nor option (1st and 2nd arg.) are applied, only 1 event is stored (3rd arg.)
      }else if(keyname == "MetaData" or keyname == "ParameterSets" or keyname == "Runs"){ // All branches from these trees are included but only 1 event is stored.
        originalTrees[keyname] = ((TTree*) rfile->Get(keyname.c_str()))->CopyTree("1","",1);
      }else if(keyname == "LuminosityBlocks"){
        originalTrees[keyname] = ((TTree*) rfile->Get(keyname.c_str()))->CopyTree("1");  // All events for this tree are stored since they are useful when comparing with lumi filtering (JSON)
      }else if( std::string(inkey->ClassName()) == "TTree"){
        std::cout << "Not copying unknown tree " << inkey->GetName() << std::endl;
      }
    }
    routfile->cd();
    for(auto tree : originalTrees){
      tree.second->Write();
    }
    rfile->Close();
    delete rfile;
  }

  std::cout << "Finished copying minimal original trees." << std::endl;
  std::cout << "------------------------------------------------------------ " << std::endl;
}

std::unordered_map<CUTS, std::vector<int>*, EnumHash> Analyzer::getArray() {
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new std::vector<int>();
  }
  return rmap;
}



void Analyzer::create_fillInfo() {

  fillInfo["FillLeadingJet"] = new FillVals(CUTS::eSusyCom, FILLER::Dipart, _Jet, _Jet);
  fillInfo["FillGen"] =        new FillVals(CUTS::eGen, FILLER::Single, _Gen);
  fillInfo["FillTau1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTau2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuon1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuon2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectron1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectron2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);

  fillInfo["FillJet1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJet2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);
  fillInfo["FillBJet"] =       new FillVals(CUTS::eRBJet, FILLER::Single, _Jet);
  fillInfo["FillCentralJet"] = new FillVals(CUTS::eRCenJet, FILLER::Single, _Jet);
  fillInfo["FillWJet"] =       new FillVals(CUTS::eRWjet, FILLER::Single, _FatJet);

  fillInfo["FillDiElectron"] = new FillVals(CUTS::eDiElec, FILLER::Dipart, _Electron, _Electron);
  fillInfo["FillDiMuon"] =     new FillVals(CUTS::eDiMuon, FILLER::Dipart, _Muon, _Muon);
  fillInfo["FillDiTau"] =      new FillVals(CUTS::eDiTau, FILLER::Dipart, _Tau, _Tau);
  fillInfo["FillMetCuts"] =    new FillVals();
  fillInfo["FillDiJet"] =      new FillVals(CUTS::eDiJet, FILLER::Dipart, _Jet, _Jet);

  fillInfo["FillMuon1Tau1"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon1Tau2"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau1"] =       new FillVals(CUTS::eMuon2Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau2"] =       new FillVals(CUTS::eMuon2Tau2, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillElectron1Tau1"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron1Tau2"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau1"] =   new FillVals(CUTS::eElec2Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau2"] =   new FillVals(CUTS::eElec2Tau2, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillMuon1Electron1"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon1Electron2"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron1"] =  new FillVals(CUTS::eMuon2Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron2"] =  new FillVals(CUTS::eMuon2Elec2, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillElectron1Jet1"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron1Jet2"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet1"] =   new FillVals(CUTS::eElec2Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet2"] =   new FillVals(CUTS::eElec2Jet2, FILLER::Dilepjet, _Electron, _Jet);

  //////I hate this solution so much.  Its terrible
  fillInfo["FillElectron1Electron2"] =     new FillVals(CUTS::eDiElec, FILLER::Single, _Electron, _Electron);
  fillInfo["FillMuon1Muon2"] =             new FillVals(CUTS::eDiMuon, FILLER::Single, _Muon, _Muon);
  fillInfo["FillTau1Tau2"] =               new FillVals(CUTS::eDiTau, FILLER::Single, _Tau, _Tau);

  //efficiency plots
  //In principal the efficiency plots should only be used, when also the object is used, but hey nobody knows!
  fillInfo["FillTauEfficiency1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTauEfficiency2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuonEfficiency1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuonEfficiency2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectronEfficiency1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectronEfficiency2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);
  fillInfo["FillJetEfficiency1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJetEfficiency2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);



  for(auto it: *histo.get_groups()) {
    if(fillInfo[it] == nullptr) fillInfo[it] = new FillVals();
  }

  // BRENDA JUN 11
  for(auto it: *syst_histo.get_groups()) {
    if(fillInfo[it] == nullptr) fillInfo[it] = new FillVals();
  }

}

void Analyzer::setupCR(std::string var, double val) {
  std::smatch m;
  std::regex part ("^(.+)_(.+)$");
  if(std::regex_match(var, m, part)) {
    std::string name = m[1];
    std::string cut = "Fill" + name;
    if(fillInfo.find(cut) == fillInfo.end()) {
      std::cout << cut << " not found, put into fillInfo" << std::endl;
      exit(1);
    }
    std::cout << cut << " " << m[2] << " " << val << " " << name << std::endl;
    testVec.push_back(new CRTester(fillInfo.at(cut), m[2], val, name));
  } else {
    std::cout << "Could not process line: " << var << std::endl;
    exit(1);
  }

}


////destructor
Analyzer::~Analyzer() {
  clear_values();
  delete BOOM;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  delete _FatJet;

  if(!isData){
    delete _Gen;
    delete _GenHadTau;
    delete _GenJet;
    delete _Photon;
  }

  for(auto fpair: fillInfo) {
    delete fpair.second;
    fpair.second=nullptr;
  }

  for(auto e: Enum<CUTS>()) {
    delete goodParts[e];
    goodParts[e]=nullptr;
  }
  for(auto it: testVec){
    delete it;
    it=nullptr;
  }

}


///resets values so analysis can start
void Analyzer::clear_values() {

  for(auto e: Enum<CUTS>()) {
    goodParts[e]->clear();
    goodParts[e]->shrink_to_fit();
  }
  //faster!!
  for(auto &it: syst_parts) {
    if (it.size() == 0) continue;
    for(auto e: Enum<CUTS>()) {
      it[e]->clear();
      it[e]->shrink_to_fit();
    }
  }
  if(infoFile!=BOOM->GetFile()){
    std::cout<<"New file!"<<std::endl;
    infoFile=BOOM->GetFile();
  }

  leadIndex=-1;
  maxCut = 0;
}

// New function: this sets up parameters that only need to be called once per event.
void Analyzer::setupEventGeneral(int nevent){

  // This function is an intermediate step called from preprocess that will set up all the variables that are common to the event and not particle specific.
  // We want to set those branches first here and then call BOOM->GetEntry(nevent) so that the variables change properly for each event.

  // For MC samples, set number of true pileup interactions, gen-HT and gen-weights.
  if(!isData){

    SetBranch("Pileup_nTrueInt",nTruePU);
    SetBranch("genWeight",gen_weight);

    if (BOOM->FindBranch("L1PreFiringWeight_Nom") != 0){
      SetBranch("L1PreFiringWeight_Nom", l1prefiringwgt);

      if(distats["Systematics"].bfind("useSystematics")){
        SetBranch("L1PreFiringWeight_Up", l1prefiringwgt_up);
        SetBranch("L1PreFiringWeight_Dn", l1prefiringwgt_dn);
      }
    }

    if (BOOM->FindBranch("LHE_HT") != 0){
      SetBranch("LHE_HT",generatorht);
    }

    if(BOOM->FindBranch("GenMET_pt") != 0){
      SetBranch("GenMET_pt", genmet_pt);
      SetBranch("GenMET_phi", genmet_phi);
    }

  }
  // Get the number of primary vertices, applies to both data and MC
  SetBranch("PV_npvs", totalVertices);
  SetBranch("PV_npvsGood", bestVertices);

  // Get the offset energy density for jet energy corrections: https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC
  SetBranch("fixedGridRhoFastjetAll", jec_rho);

  // Finally, call get entry so all the branches assigned here are filled with the proper values for each event.
  BOOM->GetEntry(nevent);

  // Check that the sample does not have crazy values of nTruePU
  if(nTruePU < 100.0){
       // std::cout << "pileupntrueint = " << pileupntrueint << std::endl;
  }
  else{
    // std::cout << "event with abnormal pileup = " << pileupntrueint << std::endl;
    clear_values();
    return;
  }

  // Calculate the pu_weight for this event.
  // pu_weight = (!isData && CalculatePUSystematics) ? hPU[(int)(nTruePU+1)] : 1.0;
  pu_weight = (!isData && CalculatePUSystematics) ? hist_pu_wgt->GetBinContent(hist_pu_wgt->FindBin(nTruePU)) : 1.0;
  // std::cout << "ntruepu = " << nTruePU << ", pu_weight (array) = " << pu_weight << ", bin in hist_ratio = " << (hist_pu_wgt->FindBin(nTruePU)) << ", pu weight (histo) = " << (hist_pu_wgt->GetBinContent(hist_pu_wgt->FindBin(nTruePU))) << std::endl;

  // Get the trigger decision vector.
  triggerDecision = false; // Reset the decision flag for each event.

  if(trigger1BranchesList.size() > 0){
    for(std::string triggname : trigger1BranchesList){
      // std::cout << "Trigger name: " << triggname << std::endl;

      TBranch *triggerBranch = BOOM->GetBranch(triggname.c_str());
      triggerBranch->SetStatus(1);
      triggerBranch->SetAddress(&triggerDecision);

      // SetBranch(triggname.c_str(), triggerDecision);
      BOOM->GetEntry(nevent);

      // std::cout << "Decision = " << triggerDecision << std::endl;
      trigger1namedecisions.push_back(triggerDecision);
      triggerBranch->ResetAddress();
    }
  }

  if(trigger2BranchesList.size() > 0){
    for(std::string triggname : trigger2BranchesList){
      // std::cout << "Trigger name: " << triggname << std::endl;

      TBranch *triggerBranch = BOOM->GetBranch(triggname.c_str());
      triggerBranch->SetStatus(1);
      triggerBranch->SetAddress(&triggerDecision);

      // SetBranch(triggname.c_str(), triggerDecision);
      BOOM->GetEntry(nevent);

      // std::cout << "Decision = " << triggerDecision << std::endl;
      trigger2namedecisions.push_back(triggerDecision);
      triggerBranch->ResetAddress();
    }
  }

  /*
  for(size_t i=0; i < triggernamedecisions.size(); i++){
    std::cout << "Trigger decision #" << i << " = " << triggernamedecisions.at(i) << std::endl;
  }
  */

}

bool Analyzer::passGenHTFilter(float genhtvalue){

  if(genhtvalue >= distats["Run"].dmap.at("LowerGenHtCut") && genhtvalue <= distats["Run"].dmap.at("UpperGenHtCut")){
    //std::cout << "genhtvalue = " << genhtvalue << ", passed genht filter " << std::endl;
    return true;
  }
  else{
    //std::cout << "genhtvalue = " << genhtvalue << ", failed genht filter " << std::endl;
    return false;
  }

}

bool Analyzer::passGenMassFilterZ(float mass_lowbound, float mass_upbound){

  if(!isZsample) return true;

   std::vector<uint> genleptonindices;
   int genpart_id = 0, genmotherpart_idx = 0, genmotherpart_id = 0;

   // std::cout << "---------------------------" << std::endl;
   // Loop over all generator level particles to look for the leptons that come from a Z:
   for(size_t idx = 0; idx < _Gen->size(); idx++) {
     // Get the particle PDG ID
     genpart_id = abs(_Gen->pdg_id[idx]);
     // Find the index of the mother particle
     genmotherpart_idx = _Gen->genPartIdxMother[idx];

     // Find the id of the mother particle using the index we just retrieved
     genmotherpart_id = _Gen->pdg_id[genmotherpart_idx];

     // std::cout << "part_idx = " << idx << ", genpart_id = " << genpart_id << ", status = " << _Gen->status[idx] << ", mother part_idx = " << genmotherpart_idx << ", mother genpart_id = " << genmotherpart_id << std::endl;

     // Only select those particles that are electrons (11), muons (13) or taus (15)
     if(! ( ( (genpart_id == 11 || genpart_id == 13) && _Gen->status[idx] == 1) || (genpart_id == 15 && _Gen->status[idx] == 2) ) )  continue;
     // std::cout << "This is a lepton" << std::endl;

     // Check that the mother particle is a Z boson (23)
     if(genmotherpart_id != 23) continue;
     // std::cout << "The mother is a Z boson" << std::endl;
     // Add the 4-momentum of this particle to the vector used later for dilepton mass calculation.
     genleptonindices.push_back(idx);

   }

   // Check that this vector only contains the information from two leptons.
   if(genleptonindices.size() != 2){
     // std::cout << "Not enough or too much leptons coming from the Z, returning false" << std::endl;
     return false;
   }
   // Check that they have the same mother (same mother index)
   else if(_Gen->genPartIdxMother[genleptonindices.at(0)] != _Gen->genPartIdxMother[genleptonindices.at(1)]){
     // std::cout << "The mothers of these two leptons are different, returning false" << std::endl;
     return false;
   }

   // std::cout << "We got a real Z -> ll event" << std::endl;
   // Get the dilepton mass using the indices from the leptons stored
   TLorentzVector lep1 = _Gen->p4(genleptonindices.at(0));
   TLorentzVector lep2 = _Gen->p4(genleptonindices.at(1));

   gendilepmass = (lep1 + lep2).M();

   if( gendilepmass >= mass_lowbound && gendilepmass <= mass_upbound){
     // std::cout << "Dilepton mass = " << (lep1 + lep2).M() << " between " << mass_lowbound << " and " << mass_upbound << " Passed the genMassfilter! " << std::endl;
     return true;
   }
   else{
     // std::cout << "Dilepton mass = " << (lep1 + lep2).M() << " not between " << mass_lowbound << " and " << mass_upbound << " Failed the genMassfilter! " << std::endl;
     return false;
   }

 }


bool Analyzer::checkGoodRunsAndLumis(int event){

    UInt_t run_num;  //NEW:  create run_num to store the run number.
    UInt_t luminosityBlock_num;  //NEW:  create luminosityBlock_num to store the number of the luminosity block.

    SetBranch("run",run_num);  //NEW:  define the branch for the run number.
    SetBranch("luminosityBlock",luminosityBlock_num);  //NEW:  define the branch for the luminosity block.

    BOOM->GetEntry(event);  //NEW:  get event.

    int key = run_num;
    int element = luminosityBlock_num;

    auto search = jsonlinedict.find(key);  //NEW:  see if the run number is in the json dictionary.

    if(search != jsonlinedict.end()){  //NEW:  this means that the run is in there.

		std::vector<int> lumivector;  //NEW:  going to make a container to store the lumis for the run.

		for (auto itr = jsonlinedict.begin(); itr != jsonlinedict.end(); itr++){ //NEW:  go through all the pairs of run, lumibound in the dictionary.
			if (itr->first != key) continue; //{ //NEW:  look for the run number.
			lumivector.push_back(itr->second);  //NEW:  grab all of the lumibounds corresponding to the run number.
		}

		// NEW:  going to go through pairs.  good lumi sections defined by [lumibound1, lumibound2], [lumibound3, lumibound4], etc.
		// This is why we have to step through the check in twos.
		for(size_t lowbound = 0; lowbound < lumivector.size(); lowbound = lowbound+2){
			int upbound = lowbound+1; //NEW:  checking bounds that are side-by-side in the vector.
			if(!(element >= lumivector[lowbound] && element <= lumivector[upbound])) continue; //NEW:  if the lumisection is not within the bounds continue, check all lumi pairs.
			// if the lumisection falls within one of the lumipairs, you will make it here. Then, return true
			return true;
		}
    }
    else{

	   return false;
    }

    return false;

}

bool Analyzer::passHEMveto2018(){

  bool hasJetinHEM = false;
  // Loop over all jets in the reco Jet collection before applying any selections.
  for(size_t index = 0; index < _Jet->size(); index++){
    // Get the 4-momentum of the jet:
    TLorentzVector jetP4 = _Jet->p4(index);

    //std::cout << "Jet #" << index << ": pt = " << jetP4.Pt() << ", eta = " << jetP4.Eta() << ", phi = " << jetP4.Phi() << std::endl;
    // Check if jet is in the HEM region
    if(jetP4.Pt() >= 30.0 && jetP4.Eta() >= -3.0 && jetP4.Eta() <= -1.3 && jetP4.Phi() >= -1.57 && jetP4.Phi() <= -0.87){
      //std::cout << "....... This jet is in the HEM region" << std::endl;
      hasJetinHEM = true;
      break;
    }
  }

  if(hasJetinHEM == true) return false;

  return true;

}

bool Analyzer::skimSignalMC(int event){
  // Check if we should use this function at all... to be called only for signal samples

  if(!isSignalMC) return true;

  // Construct the name of the branch:
  std::string signalBranchName = ("GenModel_"+inputSignalModel+"_"+inputSignalMassParam).c_str();

  // std::cout << "Name of the branch: " << signalBranchName << std::endl;
  bool isInputSignal = false;

  TBranch *signalBranch = BOOM->GetBranch(signalBranchName.c_str());
  signalBranch->SetStatus(1);
  signalBranch->SetAddress(&isInputSignal);

  BOOM->GetEntry(event);

  finalInputSignal = isInputSignal;

  signalBranch->ResetAddress();
  // std::cout << "Input signal: " << signalBranchName << ", isInputSignal? " << finalInputSignal << std::endl;
  return finalInputSignal;
}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event, std::string year){ // This function no longer needs to get the JSON dictionary as input.

  int test= BOOM->GetEntry(event);
  // std::cout << "Event test: " << test << std::endl;
  if(test<0){
    std::cout << "Could not read the event from the following file: "<<BOOM->GetFile()->GetNewUrl().Data() << std::endl;
  }

  runyear = year;
  for(Particle* ipart: allParticles){
    ipart->init();
  }
  _MET->init();


  active_part = &goodParts;
  additionalEEnoiseJets.clear();
  additionalEEnoiseJets.shrink_to_fit();

  // Call the new function setupEventGeneral: this will set generatorht, pu weight and genweight
  setupEventGeneral(event);

  if(!isData){ // Do everything that corresponds only to MC

    // Initialize the lists of generator-level particles.
    _Gen->setOrigReco();
    _GenHadTau->setOrigReco();
    _GenJet->setOrigReco();

    getGoodGen(_Gen->pstats["Gen"]);
    getGoodGenHadronicTaus(_GenHadTau->pstats["Gen"]);
    getGoodGenJets(_GenJet->pstats["Gen"]);
    getGoodGenBJets(_GenJet->pstats["Gen"]);
    getGoodGenHadronicTauNeutrinos(_Gen->pstats["Gen"]);
    // getGoodGenBJet(); //01.16.19

    //--- filtering inclusive HT-binned samples: must be done after setupEventGeneral --- //

    if(distats["Run"].bfind("DiscrByGenHT")){
      //std::cout << "generatorht = " << generatorht << std::endl;
      if(passGenHTFilter(generatorht) == false){
        clear_values();
        return;
      }
    }

    //--- filtering inclusive HT-binned samples: must be done after  _Gen->setOrigReco() --- //
     if(distats["Run"].bfind("DiscrByGenDileptonMass")){
       if(passGenMassFilterZ(distats["Run"].pmap.at("GenDilepMassRange").first, distats["Run"].pmap.at("GenDilepMassRange").second) == false){
         clear_values();
         return;
       }
     }

     // -- For signal samples -- //
     if(isSignalMC){
      if(skimSignalMC(event) == false){
        clear_values();
        return;
      }
     }

  }
  else if(isData){

  	// If you want to filter data by good run and lumi sections, turn on this option in Run_info.in
  	// SingleMuon data sets need to be filtered. SingleElectron does not (?).
  	if(distats["Run"].bfind("FilterDataByGoldenJSON")){
	    if(checkGoodRunsAndLumis(event) == false){
	    	clear_values();
	    	return;
	    }
  	}
  }

  // Call the new function passMetFilters
  // Apply here the MET filters, in case the option is turned on. It applies to both data and MC
  applymetfilters = distats["Run"].bfind("ApplyMetFilters");
  if(applymetfilters){
  	passedmetfilters = passMetFilters(year, event);
  	if(!passedmetfilters){
  		clear_values();
  		return;
  	}
  }

  // Apply HEM veto for 2018 if the flag is on.
  bool checkHEM = distats["Run"].bfind("ApplyHEMVeto2018");
  if(checkHEM == 1 || checkHEM == true){
    if(passHEMveto2018() == false){
      clear_values();
      return;
    }
  }

  // std::cout << "------------" << std::endl;

  // ------- Number of primary vertices requirement -------- //
  active_part->at(CUTS::eRVertex)->resize(bestVertices);

  // ---------------- Trigger requirement ------------------ //
  TriggerCuts(CUTS::eRTrig1);
  TriggerCuts(CUTS::eRTrig2);

  // ---------------- Additional EE noise veto -------------- //
  passAdditionalEENoiseVeto(CUTS::eREENoiseVeto, _Jet->pstats["Jet1"], year);
  

  for(size_t i=0; i < syst_names.size(); i++) {
  	std::string systname = syst_names.at(i);

    //////Smearing
    smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"], distats["Electron_systematics"], i);
    smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"], distats["Muon_systematics"], i);
    // smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"], distats["Tau_systematics"], i);
    smearTaus(*_Tau, _Tau->pstats["Smear"], distats["Tau_systematics"], i);
    //std::cout << "(Before smearing) systematic = " << systname << std::endl;
    //for(size_t i=0; i < _Jet->size(); i++){
    //  std::cout << "Jet #" << i << ": pt = " << _Jet->pt(i) << ", eta = " << _Jet->eta(i) << ", phi = " << _Jet->phi(i) << ", mass = " << _Jet->mass(i) << std::endl;
    //}

    applyJetEnergyCorrections(*_Jet,CUTS::eGJet,_Jet->pstats["Smear"], year, i);
    applyJetEnergyCorrections(*_FatJet,CUTS::eGJet,_FatJet->pstats["Smear"], year, i);
    updateMet(i);

  }

  for(size_t i=0; i < syst_names.size(); i++) {
    std::string systname = syst_names.at(i);
    for( auto part: allParticles) part->setCurrentP(i);
    _MET->setCurrentP(i);
    //std::cout << "(After smearing) systematic = " << systname << std::endl;
    //for(size_t i=0; i < _Jet->size(); i++){
    //  std::cout << "Jet #" << i << ": pt = " << _Jet->pt(i) << ", eta = " << _Jet->eta(i) << ", phi = " << _Jet->phi(i) << ", mass = " << _Jet->mass(i) << std::endl;
    //}

    getGoodParticles(i);
    selectMet(i);

    active_part=&syst_parts.at(i);
  }

  // std::cout << typeid(*active_part->at(CUTS::eRJet1)).name() << std::endl;

  active_part = &goodParts;

  if( event < 10 || ( event < 100 && event % 10 == 0 ) ||
  ( event < 1000 && event % 100 == 0 ) ||
  ( event < 10000 && event % 1000 == 0 ) ||
  ( event >= 10000 && event % 10000 == 0 ) ) {
    std::cout << std::setprecision(2)<<event << " Events analyzed "<< static_cast<double>(event)/nentries*100. <<"% done"<<std::endl;
    std::cout << std::setprecision(5);
  }
}


void Analyzer::getGoodParticles(int syst){

  std::string systname=syst_names.at(syst);
  if(syst == 0) active_part = &goodParts;
  else active_part=&syst_parts.at(syst);
    //    syst=syst_names[syst];



  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"],syst);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGHadTau, _Tau->pstats["Tau1"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGHadTau, _Tau->pstats["Tau2"],syst);
  getGoodRecoBJets(CUTS::eRBJet, CUTS::eGBJet, _Jet->pstats["BJet"],syst); //01.16.19
  getAdditionalBadEENoiseRecoJets(_Jet->pstats["Jet1"]);
  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"],syst);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"],syst);
  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"],syst);
  getGoodRecoLeadJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"],syst);
  getGoodRecoLeadJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"],syst);

  getGoodRecoFatJets(CUTS::eRWjet, _FatJet->pstats["Wjet"],syst);

  ///VBF Susy cut on leadin jets
  VBFTopologyCut(distats["VBFSUSY"],syst);

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1,CUTS::eRTau1, CUTS::eElec1Tau1, distats["Electron1Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau1, CUTS::eElec2Tau1, distats["Electron2Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1, CUTS::eRTau2, CUTS::eElec1Tau2, distats["Electron1Tau2"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau2, CUTS::eElec2Tau2, distats["Electron2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau1, CUTS::eMuon1Tau1, distats["Muon1Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau2, CUTS::eMuon1Tau2, distats["Muon1Tau2"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau1, CUTS::eMuon2Tau1, distats["Muon2Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau2, CUTS::eMuon2Tau2, distats["Muon2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec1, CUTS::eMuon1Elec1, distats["Muon1Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec2, CUTS::eMuon1Elec2, distats["Muon1Electron2"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec1, CUTS::eMuon2Elec1, distats["Muon2Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec2, CUTS::eMuon2Elec2, distats["Muon2Electron2"],syst);

  ////Dilepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"],syst);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"],syst);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"],syst);

  //
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec1, CUTS::eRJet1, CUTS::eElec1Jet1, distats["Electron1Jet1"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec1, CUTS::eRJet2, CUTS::eElec1Jet2, distats["Electron1Jet2"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec2, CUTS::eRJet1, CUTS::eElec2Jet1, distats["Electron2Jet1"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec2, CUTS::eRJet2, CUTS::eElec2Jet2, distats["Electron2Jet2"],syst);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"],syst);

}


void Analyzer::fill_efficiency() {
  //cut efficiency
  const std::vector<CUTS> goodGenLep={CUTS::eGElec,CUTS::eGMuon,CUTS::eGTau};
  //just the lepton 1 for now
  const std::vector<CUTS> goodRecoLep={CUTS::eRElec1,CUTS::eRMuon1,CUTS::eRTau1};



  for(size_t igen=0;igen<goodGenLep.size();igen++){
    Particle* part =particleCutMap.at(goodGenLep[igen]);
    CUTS cut=goodRecoLep[igen];
    std::smatch mGen;
    std::string tmps=part->getName();
    std::regex_match(tmps, mGen, genName_regex);
    //loop over all gen leptons
    for(int iigen : *active_part->at(goodGenLep[igen])){


      int foundReco=-1;
      for(size_t ireco=0; ireco<part->size(); ireco++){
        if(part->p4(ireco).DeltaR(_Gen->p4(iigen))<0.3){
          foundReco=ireco;
        }
      }
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),foundReco>=0,0);
      if(foundReco>=0){
        bool id_particle= (find(active_part->at(cut)->begin(),active_part->at(cut)->end(),foundReco)!=active_part->at(cut)->end());
        histo.addEffiency("eff_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),id_particle,0);
      }
    }
  }

}


////Reads cuts from Cuts.in file and see if the event has enough particles
bool Analyzer::fillCuts(bool fillCounter) {
  const std::unordered_map<std::string,std::pair<int,int> >* cut_info = histo.get_cuts();
  const std::vector<std::string>* cut_order = histo.get_cutorder();

  bool prevTrue = true;

  maxCut=0;

  for(size_t i = 0; i < cut_order->size(); i++) {
    std::string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != std::string::npos){
      maxCut += 1;
      continue;
    }
    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = active_part->at(cut_num.at(cut))->size();

    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num.at(cut) == CUTS::eR1stJet || cut_num.at(cut) == CUTS::eR2ndJet) && active_part->at(cut_num.at(cut))->at(0) == -1 ) {
        prevTrue = false;
        continue;  ////dirty dirty hack
      }
      if(fillCounter && crbins == 1) {
        cuts_per[i]++;
        cuts_cumul[i] += (prevTrue) ? 1 : 0;
        maxCut += (prevTrue) ? 1 : 0;
      }else{
        maxCut += (prevTrue) ? 1 : 0;
      }
    }else {
      prevTrue = false;
    }
  }

  if(crbins != 1) {
    if(!prevTrue) {
      maxCut = -1;
      return prevTrue;
    }

    int factor = crbins;
    for(auto tester: testVec) {
      factor /= 2;
      /////get variable value from maper.first.
      if(tester->test(this)) { ///pass cut
        maxCut += factor;
      }
    }
    if(isData && blinded && maxCut == SignalRegion) return false;
    cuts_per[maxCut]++;
  }

  return prevTrue;
}



///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  std::vector<std::string> cut_order;
  if(crbins > 1) cut_order = *(histo.get_folders());
  else cut_order = *(histo.get_cutorder());
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  double run_time_real=elapsed_seconds.count();


  std::cout.setf(std::ios::floatfield,std::ios::fixed);
  std::cout<<std::setprecision(3);
  std::cout << "\n";
  std::cout << "Selection Efficiency " << "\n";
  std::cout << "Total events: " << nentries << "\n";
  std::cout << "\n";
  std::cout << "Run Time (real): " <<run_time_real <<" s\n";
  std::cout << "Time per 1k Events (real): " << run_time_real/(nentries/1000) <<" s\n";
  std::cout << "Events/s: " << static_cast<double>(nentries)/(run_time_real) <<" 1/s (real) \n";
  std::cout << "                        Name                  Indiv.";
  if(crbins == 1) std::cout << "            Cumulative";
  std::cout << std::endl << "---------------------------------------------------------------------------\n";
  for(size_t i = 0; i < cut_order.size(); i++) {
    std::cout << std::setw(28) << cut_order.at(i) << "    ";
    if(isData && cut_order.at(i).find("Gen") != std::string::npos) std::cout << "Skipped" << std::endl;
    else if(crbins != 1 && blinded && i == (size_t)SignalRegion) std::cout << "Blinded Signal Region" << std::endl;
    else {
      std::cout << std::setw(10) << cuts_per.at(i) << "  ( " << std::setw(5) << ((float)cuts_per.at(i)) / nentries << ") ";
      if(crbins == 1) std::cout << std::setw(12) << cuts_cumul.at(i) << "  ( " << std::setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") ";

      std::cout << std::endl;
    }
  }
  std::cout <<std::setprecision(5);
  std::cout << "---------------------------------------------------------------------------\n";

  //write all the histograms
  //attention this is not the fill_histogram method from the Analyser
  histo.fill_histogram(routfile);
  if(doSystematics)
    syst_histo.fill_histogram(routfile);

  routfile->Write();
  routfile->Close();

}

/////////////PRIVATE FUNCTIONS////////////////


bool Analyzer::select_mc_background(){
  //will return true if Z* mass is smaller than 200GeV
  if(_Gen == nullptr){
    return true;
  }
  if(gen_selection["DY_noMass_gt_200"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<200;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //will return true if Z* mass is smaller than 200GeV
  if(gen_selection["DY_noMass_gt_100"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<100;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //cout<<"Something is rotten in the state of Denmark."<<std::endl;
  //cout<<"could not find gen selection particle"<<std::endl;
  return true;
}

// --- Function that sets up all the SFs recommended by the Tau POG for Run II legacy --- //
void Analyzer::setupTauIDSFsInfo(std::string tauidalgoname, std::string year, bool applyDMsfs, bool applyEmbedding){

  static std::map<std::string, std::string> tauidyearmap = {
    {"2016", "2016Legacy"},
    {"2017", "2017ReReco"},
    {"2018", "2018ReReco"}
  };

  tauidyear = tauidyearmap[year];

  // Read the corresponding ID algorithm:
  if(tauidalgoname.find("DeepTau") != std::string::npos){
    tauid_algo = "DeepTau2017v2p1VSjet";
    antiele_algo = "DeepTau2017v2p1VSe";
    antimu_algo = "DeepTau2017v2p1VSmu";

    tauidwpsmap = {
      {1, "VVVLoose"}, {2, "VVLoose"}, {4, "VLoose"}, {8, "Loose"},
      {16, "Medium"}, {32, "Tight"}, {64, "VTight"}, {128, "VVTight"}
    };

    antielewpsmap = {
      {1, "VVVLoose"}, {2, "VVLoose"}, {4, "VLoose"}, {8, "Loose"},
      {16, "Medium"}, {32, "Tight"}, {64, "VTight"}, {128, "VVTight"}
    };

    antimuwpsmap = {
      {1, "VLoose"}, {2, "Loose"}, {4, "Medium"}, {8, "Tight"}
    };

  }
  else if(tauidalgoname.find("MVA") != std::string::npos){
    tauid_algo = "MVAoldDM2017v2";
    antiele_algo = "antiEleMVA6";
    antimu_algo = "antiMu3";

    tauidwpsmap = {
      {1, "VVLoose"}, {2, "VLoose"}, {4, "Loose"}, {8, "Medium"},
      {16, "Tight"}, {32, "VTight"}, {64, "VVTight"}
    };

    antielewpsmap = {
      {1, "VLoose"}, {2, "Loose"}, {4, "Medium"}, {8, "Tight"}, {16, "VTight"}
    };

    antimuwpsmap = {
      {1, "Loose"}, {2, "Tight"}
    };
  }

  //std::cout << "Tau IDs: " << tauid_algo << ", " << antiele_algo << ", " << antimu_algo << std::endl;

  // ----- Tau 1 stuff ----- //
  if(!_Tau->pstats["Tau1"].bfind("FlipIsolationRequirement")){
    tauidwp = _Tau->pstats["Tau1"].dmap.at("DiscrByMinIsolation");
  }
  else{
    tauidwp = _Tau->pstats["Tau1"].dmap.at("DiscrByMaxIsolation");
  }
  antielewp = _Tau->pstats["Tau1"].dmap.at("DiscrAgainstElectron");
  antimuwp = _Tau->pstats["Tau1"].dmap.at("DiscrAgainstMuon");

  //std::cout << "Working points 1: tau id = " << tauidwp << ", antiele = " << antielewp << ", antimu = " << antimuwp << std::endl;
  //std::cout << "Working points 1: tau id = " << tauidwpsmap[tauidwp] << ", antiele = " << antielewpsmap[antielewp] << ", antimu = " << antimuwpsmap[antimuwp] << std::endl;

  // Load the modules according to the algorithm and working points specified for tau ID SFs
  if(tauidwp != 0){
    tau1idSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, tauid_algo, tauidwpsmap[tauidwp], applyDMsfs, applyEmbedding);
  }
  else{
    failtau1iso = true;
  }
  tau1id_antiEleSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, antiele_algo, antielewpsmap[antielewp]);
  tau1id_antiMuSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, antimu_algo, antimuwpsmap[antimuwp]);


  // Reset all wp variables
  tauidwp = 0;
  antielewp = 0;
  antimuwp = 0;


  // ------  Tau 2 stuff  ------ //
  if(!_Tau->pstats["Tau2"].bfind("FlipIsolationRequirement")){
    tauidwp = _Tau->pstats["Tau2"].dmap.at("DiscrByMinIsolation");
  }
  else{
    tauidwp = _Tau->pstats["Tau2"].dmap.at("DiscrByMaxIsolation");
  }
  antielewp = _Tau->pstats["Tau2"].dmap.at("DiscrAgainstElectron");
  antimuwp = _Tau->pstats["Tau2"].dmap.at("DiscrAgainstMuon");

  //std::cout << "Working points 2: tau id = " << tauidwp << ", antiele = " << antielewp << ", antimu = " << antimuwp << std::endl;

  // Load the modules according to the algorithm and working points specified for tau ID SFs
  if(tauidwp != 0){
    tau2idSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, tauid_algo, tauidwpsmap[tauidwp], applyDMsfs, applyEmbedding);
  }
  else{
    failtau2iso = true;
  }
  tau2id_antiEleSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, antiele_algo, antielewpsmap[antielewp]);
  tau2id_antiMuSFs = TauIDSFTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, antimu_algo, antimuwpsmap[antimuwp]);

}

void Analyzer::setupTauResSFsInfo(bool taufakeenergyscale){

  // Load the corresponding modules for tau energy scale and tau fake energy scale
  tauesSFs = TauESTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, tauid_algo);

  if(antiele_algo.find("DeepTau") != std::string::npos && taufakeenergyscale){
    taufesSFs = TauFESTool((PUSPACE+"TauIDSFs/data").c_str(), tauidyear, antiele_algo);
  }

}

double Analyzer::getTauIdSFs(bool getTauIDsf, bool getTauIDbyDMsfs, bool getAntiElesf, bool getAntiMusf, std::string uncertainty){

  double sf = 1.0, sf_tauid = 1.0, sf_antiele = 1.0, sf_antimu = 1.0;

  // std::cout << "Size of Tau1 collection = " << active_part->at(CUTS::eRTau1)->size() << ", Tau2 collection = " << active_part->at(CUTS::eRTau2)->size() << std::endl;
  // Loop over the Tau1 list first.
  for(auto i : *active_part->at(CUTS::eRTau1)){
    // Look for the current muon in the eRMuon2 list. If found, skip it.
    if(find(active_part->at(CUTS::eRTau2)->begin(), active_part->at(CUTS::eRTau2)->end(), i) != active_part->at(CUTS::eRTau2)->end() ) continue;

    // Get the gen-match status of the current tau
    int gen_match_status = static_cast<int>(_Tau->genPartFlav[i]);

    // std::cout << "gen match status = " << gen_match_status << std::endl;

    if(getTauIDsf && !failtau1iso){
      if(!getTauIDbyDMsfs){ // pT-dependent SFs
        sf_tauid *= tau1idSFs.getSFvsPT(_Tau->pt(i), gen_match_status, uncertainty);
      }
      else{ // DM-dependent SFs
        sf_tauid *= tau1idSFs.getSFvsDM(_Tau->pt(i), _Tau->decayModeInt[i], gen_match_status, uncertainty);
      }
    }
    else if(getTauIDsf && failtau1iso) sf_tauid = 1.0;
    else if(getAntiElesf) sf_antiele *= tau1id_antiEleSFs.getSFvsEta(_Tau->eta(i), gen_match_status, uncertainty);
    else if(getAntiMusf) sf_antimu *= tau1id_antiMuSFs.getSFvsEta(_Tau->eta(i), gen_match_status, uncertainty);
  }

  // std::cout << "sf_tauid = " << sf_tauid << ", sf_antiele = " << sf_antiele << ", sf_antimu = " << sf_antimu << std::endl;

  // Now, we loop over the tau2 list, to apply the SFs to the uncorrected taus from the previous list.
  for(auto i : *active_part->at(CUTS::eRTau2)){
    // Get the gen-match status of the current tau
    int gen_match_status = static_cast<int>(_Tau->genPartFlav[i]);
    // std::cout << "gen match status = " << gen_match_status << std::endl;

    if(getTauIDsf && !failtau2iso){
      if(!getTauIDbyDMsfs){ // pT-dependent SFs
        sf_tauid *= tau2idSFs.getSFvsPT(_Tau->pt(i), gen_match_status, uncertainty);
      }
      else{ // DM-dependent SFs
        sf_tauid *= tau2idSFs.getSFvsDM(_Tau->pt(i), _Tau->decayModeInt[i], gen_match_status, uncertainty);
      }
    }
    else if(getTauIDsf && failtau2iso) sf_tauid = 1.0;
    else if(getAntiElesf) sf_antiele *= tau2id_antiEleSFs.getSFvsEta(_Tau->eta(i), gen_match_status, uncertainty);
    else if(getAntiMusf) sf_antimu *= tau2id_antiMuSFs.getSFvsEta(_Tau->eta(i), gen_match_status, uncertainty);
  }

  //std::cout << "sf_tauid = " << sf_tauid << ", sf_antiele = " << sf_antiele << ", sf_antimu = " << sf_antimu << std::endl;

  sf = sf_tauid * sf_antiele * sf_antimu;

  return sf;
}


// ---- Function that propagates the unclustered energy uncertainty to MET ---- //
void Analyzer::updateMet(int syst){

  std::string systname = syst_names.at(syst);

  // After jet energy corrections, we need to update all MET vectors in each systematics so they are the same as that in "orig"
  if(_MET->needSys(syst) != 0){

    TLorentzVector met_p4_nom = _MET->getNominalP();
    _MET->addP4Syst(met_p4_nom, syst);

  }
  // --- Apply MET unclustered energy uncertainty if needed (syst) ---- //
  if(!isData){
    _MET->propagateUnclEnergyUncty(systname,syst);
  }

}


///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::selectMet(int syst) {

  // Before using the TreatMuonsAsNeutrinos option, we store the MET value in a separate vector for reference purposes:
  _MET->JERCorrMet.SetPxPyPzE(_MET->px(), _MET->py(), _MET->p4().Pz(), _MET->energy());
  // std::cout << "Met before cuts: px = " << _MET->px() << ", py = " << _MET->py() << std::endl;
  if(distats["Run"].bfind("TreatMuonsAsNeutrinos") || distats["Run"].bfind("TreatOnlyOneMuonAsNeutrino") ) treatMuonsAsMet(syst);

  _MET->calculateHtAndMHt(distats["Run"], *_Jet,  syst);
  /////MET CUTS

  if(!passCutRange(_MET->pt(), distats["Run"].pmap.at("MetCut"))) return;

  if(distats["Run"].bfind("DiscrByHT") && _MET->HT() < distats["Run"].dmap.at("HtCut")) return;

  if(syst==0){
    active_part->at(CUTS::eMET)->push_back(1);
  }
  else{
    syst_parts.at(syst).at(CUTS::eMET)->push_back(1);
  }
}

bool Analyzer::passMetFilters(std::string year, int ievent){

  // Set the branches accordingly:
  // good vertices filter
  SetBranch("Flag_goodVertices", primaryvertexfilter);
  // beam halo filter
  SetBranch("Flag_globalSuperTightHalo2016Filter", beamhalofilter);
  // HBHE noise filter
  SetBranch("Flag_HBHENoiseFilter", hbhenoisefilter);
  // ECAL trigger primitives filter
  SetBranch("Flag_EcalDeadCellTriggerPrimitiveFilter", ecaltpfilter);
  // Bad PF muon filter
  SetBranch("Flag_BadPFMuonFilter", badpfmuonfilter);
  // Bad charged hadron filter - not recommended.
  // SetBranch("Flag_BadChargedCandidateFilter", badchargedhadronfilter);
  if(year.compare("2016") == 0){
    // in 2016, this filter is not recommended... therefore we set it always to true.
    ecalbadcalibrationfilter = true;
  }
  else{
    // ECAL bad calibration filter (2017 + 2018).
    // SetBranch("Flag_ecalBadCalibFilter", ecalbadcalibrationfilter);
    SetBranch("Flag_ecalBadCalibFilterV2", ecalbadcalibrationfilter);
  }

  // Call get entry so all the branches assigned here are filled with the proper values for each event.
  BOOM->GetEntry(ievent);

  // Check if the current event passed all the flags
  allmetfilters = primaryvertexfilter && beamhalofilter && hbhenoisefilter && ecaltpfilter && badpfmuonfilter && ecalbadcalibrationfilter;

  return allmetfilters;

}

void Analyzer::treatMuonsAsMet(int syst) {


  if( ! ( distats["Run"].bfind("TreatMuonsAsNeutrinos") ||  distats["Run"].bfind("TreatOnlyOneMuonAsNeutrino") ) ) return;

  //  Initialize the deltaMus before calculation
  _MET->systdeltaMEx[syst] = 0.0;
  _MET->systdeltaMEy[syst] = 0.0;

  if(! distats["Run"].bfind("TreatOnlyOneMuonAsNeutrino")){ // All muons in the event are being treated as neutrinos.

    // First, loop over the reco muon1 list
    for(auto it : *active_part->at(CUTS::eRMuon1)) {
      // Look for the current muon in the eRMuon2 list. If found, skip it.

      if(find(active_part->at(CUTS::eRMuon2)->begin(), active_part->at(CUTS::eRMuon2)->end(), it) != active_part->at(CUTS::eRMuon2)->end() ) continue;

      _MET->systdeltaMEx[syst] += _Muon->p4(it).Px();
      _MET->systdeltaMEy[syst] += _Muon->p4(it).Py();

    }

    // Next, loop over the reco muon2 list
    for(auto it : *active_part->at(CUTS::eRMuon2)){

      _MET->systdeltaMEx[syst] += _Muon->p4(it).Px();
      _MET->systdeltaMEy[syst] += _Muon->p4(it).Py();

    }
  }
  else if(distats["Run"].bfind("TreatOnlyOneMuonAsNeutrino") && (active_part->at(CUTS::eRMuon1)->size() > 0 || active_part->at(CUTS::eRMuon2)->size() > 0)){ // Only one muon in the event is treated as neutrino. Particular for single lepton final states. The muon gets randomly chosen.

    // Define a random generator
    TRandom *rndm1 = new TRandom3(0);
    TRandom *rndm2 = new TRandom3(0);

    size_t mu1size = active_part->at(CUTS::eRMuon1)->size();
    size_t mu2size = active_part->at(CUTS::eRMuon2)->size();

    // Define which list of muons will be looped over
    unsigned int rnd_num = rndm1->Integer(123456789);
    unsigned int rnd_muon = rndm2->Integer(987654321);

    if( ( (rnd_num % rnd_muon) % 2 == 0 && mu1size > 0) || ( (rnd_num % rnd_muon) % 2 == 1 && (mu2size == 0 && mu1size > 0) ) ){

      // The muon selected will be from the muon1 list. We select it randomly, using the size of this list as input.
      if(mu1size > 0){
        rnd_muon = rnd_muon % mu1size;
      }
      else if (mu2size > 0){
        rnd_muon = rnd_muon % mu2size;
      }

      int it = active_part->at(CUTS::eRMuon1)->at(rnd_muon);

      _MET->systdeltaMEx[syst] += _Muon->p4(it).Px();
      _MET->systdeltaMEy[syst] += _Muon->p4(it).Py();

    }
    else if( ((rnd_num % rnd_muon) % 2 == 1 && mu2size > 0) || ((rnd_num % rnd_muon) % 2 == 0 && (mu1size == 0 && mu2size > 0) ) ){
      // The muon selected will be from the muon2 list. We select it randomly, using the size of this list as input.

      if(mu2size > 0){
        rnd_muon = rnd_muon % mu2size;
      }
      else if (mu1size > 0){
        rnd_muon = rnd_muon % mu1size;
      }

      int it = active_part->at(CUTS::eRMuon2)->at(rnd_muon);

      _MET->systdeltaMEx[syst] += _Muon->p4(it).Px();
      _MET->systdeltaMEy[syst] += _Muon->p4(it).Py();
    }

  }

  // recalculate MET, adding the muon Px and Py accordingly
  _MET->update(syst);

}

/////Check if a given branch is not found in the file

void Analyzer::branchException(std::string branch){
  if(BOOM->FindBranch(branch.c_str()) == 0 ){
     throw "Branch not found in the current sample. Check the config files associated to this branch.";
  }
}

void Analyzer::getTriggerBranchesList(CUTS ePos, std::string trigger, bool usewildcard){

  if(ePos == CUTS::eRTrig1){
    if(usewildcard){
      for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){
        std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
        if(branch_name.find(trigger) == std::string::npos) continue;
        //std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
        trigger1BranchesList.push_back(branch_name);
      }
    } else {
      for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){
        std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
        // Look for branches that match exactly the trigger name
        if(branch_name.compare(trigger.c_str()) != 0) continue;
        // std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
        trigger1BranchesList.push_back(branch_name);
      }
    }

    if(trigger1BranchesList.size() == 0) throw "no branches matching this name were found. Check the trigger name requirement.\n";
  }

  if(ePos == CUTS::eRTrig2){
    if(usewildcard){
      for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){
        std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
        if(branch_name.find(trigger) == std::string::npos) continue;
        //std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
        trigger2BranchesList.push_back(branch_name);
      }
    } else {
      for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){
        std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
        // Look for branches that match exactly the trigger name
        if(branch_name.compare(trigger.c_str()) != 0) continue;
        // std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
        trigger2BranchesList.push_back(branch_name);
      }
    }

    if(trigger2BranchesList.size() == 0) throw "no branches matching this name were found. Check the trigger name requirement.\n";
  }
}

/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(std::string year) {

  genMaper = {
    {5, new GenFill(2, CUTS::eGJet)},     {6,  new GenFill(2, CUTS::eGTop)},
    {11, new GenFill(1, CUTS::eGElec)},   {13, new GenFill(1, CUTS::eGMuon)},
    {15, new GenFill(2, CUTS::eGTau)},    {23, new GenFill(62, CUTS::eGZ)},
    {24, new GenFill(62, CUTS::eGW)},      {25, new GenFill(2, CUTS::eGHiggs)},
    {5, new GenFill(52, CUTS::eGBJet)}
  };

  isData=true;
  if(BOOM->FindBranch("Pileup_nTrueInt")!=0){
    isData=false;
  }

  read_info(filespace + "ElectronTau_info.in");
  read_info(filespace + "MuonTau_info.in");
  read_info(filespace + "MuonElectron_info.in");
  read_info(filespace + "DiParticle_info.in");
  read_info(filespace + "ElectronJet_info.in");
  read_info(filespace + "VBFCuts_info.in");
  read_info(filespace + "Run_info.in");
  read_info(filespace + "Systematics_info.in");
  read_info(filespace + "SignalMC_info.in");

  // Call readinJSON and save this in the jsonlinedict we declared in Analyzer.h
  jsonlinedict = readinJSON(year);

  if(inputTrigger1Names.size() > 0){
    for(std::string trigger : inputTrigger1Names){
      try{
        getTriggerBranchesList(CUTS::eRTrig1, trigger, distats["Run"].bfind("UseTriggerWildcard"));
      }
      catch(const char* msg){
        std::cout << "ERROR! Trigger " << trigger << ": " << msg << std::endl;
        continue;
      }
    }

  // Check that there are no elements on the trigger list that refer to the same trigger to speed up the process.
  removeDuplicates(trigger1BranchesList);

  }

  if(inputTrigger2Names.size() > 0){
    for(std::string trigger : inputTrigger2Names){
      try{
        getTriggerBranchesList(CUTS::eRTrig2, trigger, distats["Run"].bfind("UseTriggerWildcard"));
      }
      catch(const char* msg){
        std::cout << "ERROR! Trigger " << trigger << ": " << msg << std::endl;
        continue;
      }
    }

    // Check that there are no elements on the trigger list that refer to the same trigger to speed up the process.
    removeDuplicates(trigger2BranchesList);
  }

  std::cout << " ---------------------------------------------------------------------- " << std::endl;
  if(trigger1BranchesList.size() > 0){
    std::cout << "Full list of trigger to be probed (1): " << std::endl;
    for(std::string name : trigger1BranchesList){
      std::cout << name << std::endl;
    }
  }
  if(trigger2BranchesList.size() > 0){
    std::cout << "Full list of trigger to be probed (2): " << std::endl;
    for(std::string name : trigger2BranchesList){
      std::cout << name << std::endl;
    }
  }
  std::cout << " ---------------------------------------------------------------------- " << std::endl;

  // Check if it is a signal MC sample:
  // isSignalMC = distats["SignalMC"].bfind("isSignalMC");

  // double check
  if(BOOM->FindBranch( ("GenModel_"+inputSignalModel+"_"+inputSignalMassParam).c_str()) == 0){
   isSignalMC = false;
   std::cout << "This is not a signal MC sample." << std::endl;
  }
  else if(BOOM->FindBranch( ("GenModel_"+inputSignalModel+"_"+inputSignalMassParam).c_str()) != 0){
   isSignalMC = true;
   std::cout << "This is a signal MC sample!" << std::endl;
  }
  std::cout << " ---------------------------------------------------------------------- " << std::endl;

  // Check if it is a signal MC sample:
  // isSignalMC = distats["SignalMC"].bfind("isSignalMC");

  // double check
  if(BOOM->FindBranch( ("GenModel_"+inputSignalModel+"_"+inputSignalMassParam).c_str()) == 0){
    isSignalMC = false;
    std::cout << "This is not a signal MC sample." << std::endl;
  }
  else if(BOOM->FindBranch( ("GenModel_"+inputSignalModel+"_"+inputSignalMassParam).c_str()) != 0){
    isSignalMC = true;
    std::cout << "This is a signal MC sample!" << std::endl;
  }
  std::cout << " ---------------------------------------------------------------------- " << std::endl;
}


void Analyzer::initializeMCSelection(std::vector<std::string> infiles) {
    // check if we need to make gen level cuts to cross clean the samples:

  isVSample = false;
  if(infiles[0].find("DY") != std::string::npos){
    isVSample = true;
    if(infiles[0].find("DYJetsToLL_M-50_HT-") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
      //gen_selection["DY_noMass_gt_200"]=true;
    //get the DY1Jet DY2Jet ...
    }else if(infiles[0].find("JetsToLL_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }

    if(infiles[0].find("DYJetsToLL_M-50_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }
  }else{
    //set it to false!!
    gen_selection["DY_noMass_gt_200"]=false;
    gen_selection["DY_noMass_gt_100"]=false;
  }

  if(infiles[0].find("WJets") != std::string::npos){
    isVSample = true;
  }
}


///parsing method that gets info on diparts and basic run info
//put in std::map called "distats"
void Analyzer::read_info(std::string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  std::ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }


  std::string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    std::vector<std::string> stemp;

    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      std::cout << "error in " << filename << "; no groups specified for data" << std::endl;
      exit(1);
    } else if(stemp.size() == 2) {
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(group.compare("Control_Region") !=0 ){
        if(stemp[1] == "1" || stemp[1] == "true"){
          distats[group].bset.push_back(stemp[0]);
          if(stemp[1] == "1" ){
            distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
          }
        }
        else if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
      }else{
        if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
      }
      if(stemp.at(0).find("Trigger1") != std::string::npos) {
        for(auto trigger : stemp){
          if(trigger.find("Trigger")== std::string::npos and "="!=trigger ){
            inputTrigger1Names.push_back(trigger);
          }
        }
        continue;
      }
      if(stemp.at(0).find("Trigger2") != std::string::npos) {
        for(auto trigger : stemp){
          if(trigger.find("Trigger")== std::string::npos and "="!=trigger ){
            inputTrigger2Names.push_back(trigger);
          }
        }
        continue;
      }
      if(stemp.at(0).find("ModelTag") != std::string::npos){
        for(auto model: stemp){
          if(model.find("ModelTag") == std::string::npos and "=" != model){
            inputSignalModel = model;
          }
        }
        continue;
      }
      if(stemp.at(0).find("MassParameters") != std::string::npos){
        for(auto massparam: stemp){
          if(massparam.find("MassParameters") == std::string::npos and "=" != massparam){
            inputSignalMassParam = massparam;
          }
        }
        continue;
      }
    } else if(stemp.size() == 3 and stemp.at(0).find("Trigger") == std::string::npos){
      distats[group].pmap[stemp[0]] = std::make_pair(std::stod(stemp[1]), std::stod(stemp[2]));
    } else{

      if(stemp.at(0).find("Trigger") != std::string::npos) {
        for(auto trigger : stemp){
          if(trigger.find("Trigger") == std::string::npos and "=" != trigger ){
            inputTrigger1Names.push_back(trigger);
          }
        }
        continue;
      }
    }
  }
  info_file.close();
}


// This code works pretty much (at least in my tests), but dagnabit, its ugly.  They all can't be winners, at least now...
void Analyzer::setCutNeeds() {


  for(auto it: *histo.get_groups()) {
    if(fillInfo[it]->type == FILLER::None) continue;
    neededCuts.loadCuts(fillInfo[it]->ePos);
  }
  for(auto it : *histo.get_cutorder()) {
    try{
      neededCuts.loadCuts(cut_num.at(it));
    }catch(...){
      std::cout<<"The following cut is strange: "<<it<<std::endl;
      exit(2);
    }
  }
  for(auto it: testVec) {
    neededCuts.loadCuts(it->info->ePos);
  }

  if(!isData and (distats["Run"].bfind("ApplyISRZBoostSF") || distats["Run"].bfind("ApplySUSYZBoostSF") || distats["Run"].bfind("ApplyVBFSusyZBoostSF")) and isVSample){
    neededCuts.loadCuts(CUTS::eGen);
    neededCuts.loadCuts(CUTS::eGZ);
    neededCuts.loadCuts(CUTS::eGW);
  }

  neededCuts.loadCuts(CUTS::eGHadTau);
  neededCuts.loadCuts(CUTS::eGBJet); //01.16.19

  neededCuts.loadCuts(_Jet->findExtraCuts());
  if(doSystematics) {
    neededCuts.loadCuts(CUTS::eGen);
  }

  for(auto it: jetCuts) {
    if(!neededCuts.isPresent(it)) continue;
    neededCuts.loadCuts(_Jet->overlapCuts(it));
  }

  if(neededCuts.isPresent(CUTS::eRWjet)) {
    neededCuts.loadCuts(_FatJet->findExtraCuts());
    neededCuts.loadCuts(_FatJet->overlapCuts(CUTS::eRWjet));
  } else {
    std::cout<<"WJets not needed. They will be deactivated!"<<std::endl;
    _FatJet->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRTau1) || neededCuts.isPresent(CUTS::eRTau2) ) {
    neededCuts.loadCuts(_Tau->findExtraCuts());
  } else {
    std::cout<<"Taus not needed. They will be deactivated!"<<std::endl;
    _Tau->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRElec1) || neededCuts.isPresent(CUTS::eRElec2) ) {
    neededCuts.loadCuts(_Electron->findExtraCuts());
  } else {
    std::cout<<"Electrons not needed. They will be deactivated!"<<std::endl;
    _Electron->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRMuon1) || neededCuts.isPresent(CUTS::eRMuon2) ) {
    neededCuts.loadCuts(_Muon->findExtraCuts());
  } else {
    std::cout<<"Muons not needed. They will be deactivated!"<<std::endl;
    _Muon->unBranch();
  }

  if( !neededCuts.isPresent(CUTS::eGen) and !isData) {
    std::cout<<"Gen not needed. They will be deactivated!"<<std::endl;
    _Gen->unBranch();

  }

  std::cout << "Cuts being filled: " << std::endl;
  for(auto cut : neededCuts.getCuts()) {
    std::cout << enumNames.at(static_cast<CUTS>(cut)) << "   ";
  }
  std::cout << std::endl;
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz std::vectors
//of the data into the std::vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lep, CUTS eGenPos, const PartStats& stats, const PartStats& syst_stats, int syst) {
  if( isData) {
    lep.setOrigReco();
    return;
  }

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) return;

  if(systname=="orig" && !stats.bfind("SmearTheParticle")){
    lep.setOrigReco();
  } else {
    systematics.loadScaleRes(stats, syst_stats, systname);
    for(size_t i = 0; i < lep.size(); i++) {
      TLorentzVector lepReco = lep.RecoP4(i);
      TLorentzVector genVec =  matchLeptonToGen(lepReco, lep.pstats["Smear"],eGenPos);
      systematics.shiftLepton(lep, lepReco, genVec, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst);
    }
  }
}

void Analyzer::smearTaus(Lepton& lep, const PartStats& stats, const PartStats& syst_stats, int syst){
  if(isData) {
    lep.setOrigReco();
    return;
  }

  if(!lep.needSyst(syst)) return;

  std::string systname = syst_names.at(syst);

  if(!stats.bfind("SmearTheParticle") && systname=="orig"){
    lep.setOrigReco();
    return;
  }

  for(size_t i = 0; i < lep.size(); i++){

    TLorentzVector tauRecoP4 = lep.RecoP4(i);

    float tes_SF = 1.0, tfes_SF = 1.0;
    double tau_pt_tesShift = tauRecoP4.Pt() * tes_SF, tau_mass_tesShift = tauRecoP4.M() * tes_SF;

    int taugenmatch = static_cast<int>(_Tau->genPartFlav[i]);

    if(stats.bfind("SmearTheParticle")){
      if(systname == "orig"){
        tes_SF = tauesSFs.getTES(tauRecoP4.Pt(), _Tau->decayModeInt[i], taugenmatch, "");
      }
      else{
        if(systname == "Tau_Scale_Up"){
          tes_SF = tauesSFs.getTES(tauRecoP4.Pt(), _Tau->decayModeInt[i], taugenmatch, "Up");
        }
        else if(systname == "Tau_Scale_Down"){
          tes_SF = tauesSFs.getTES(tauRecoP4.Pt(), _Tau->decayModeInt[i], taugenmatch, "Down");
        }
      }
    }

    if(stats.bfind("ApplyETauFakeRateESSF") && tauid_algo.find("DeepTau") != std::string::npos){
      if(systname == "orig"){
        tfes_SF = taufesSFs.getFES(tauRecoP4.Eta(), _Tau->decayModeInt[i], taugenmatch, "");
      }
      else{
        if(systname == "ETauFake_Scale_Up"){
          tfes_SF = taufesSFs.getFES(tauRecoP4.Eta(), _Tau->decayModeInt[i], taugenmatch, "Up");
        }
        else if(systname == "ETauFake_Scale_Down"){
          tfes_SF = taufesSFs.getFES(tauRecoP4.Eta(), _Tau->decayModeInt[i], taugenmatch, "Down");
        }
      }
    }

    tau_pt_tesShift = tauRecoP4.Pt() * tes_SF * tfes_SF;
    tau_mass_tesShift = tauRecoP4.M() * tes_SF * tfes_SF;

    // Shift the 4-momentum of the tau accordingly
    systematics.shiftParticle(lep, tauRecoP4, tau_pt_tesShift, tau_mass_tesShift, systname, syst);
  }

}

void Analyzer::setupJetCorrections(std::string year, std::string outputfilename){

   // ------------------------ NEW: Jet Energy Scale and Resolution corrections initialization ------------------- //
   static std::map<std::string, std::string> jecTagsMC = {
     {"2016" , "Summer16_07Aug2017_V11_MC"},
     {"2017" , "Fall17_17Nov2017_V32_MC"},
     {"2018" , "Autumn18_V19_MC"}
   };

   static std::map<std::string, std::string> jecTagsFastSim = {
     {"2016" , "Summer16_FastSimV1_MC"},
     {"2017" , "Fall17_FastSimV1_MC"},
     {"2018" , "Autumn18_FastSimV1_MC"}
   };

   static std::map<std::string, std::string> archiveTagsDATA = {
     {"2016" , "Summer16_07Aug2017_V11_DATA"},
     {"2017" , "Fall17_17Nov2017_V32_DATA"},
     {"2018" , "Autumn18_V19_DATA"}
   };

   static std::map<std::string, std::string> jecTagsDATA = {
     {"2016B" , "Summer16_07Aug2017BCD_V11_DATA"},
     {"2016C" , "Summer16_07Aug2017BCD_V11_DATA"},
     {"2016D" , "Summer16_07Aug2017BCD_V11_DATA"},
     {"2016E" , "Summer16_07Aug2017EF_V11_DATA"},
     {"2016F" , "Summer16_07Aug2017EF_V11_DATA"},
     {"2016G" , "Summer16_07Aug2017GH_V11_DATA"},
     {"2016H" , "Summer16_07Aug2017GH_V11_DATA"},
     {"2017B" , "Fall17_17Nov2017B_V32_DATA"},
     {"2017C" , "Fall17_17Nov2017C_V32_DATA"},
     {"2017D" , "Fall17_17Nov2017DE_V32_DATA"},
     {"2017E" , "Fall17_17Nov2017DE_V32_DATA"},
     {"2017F" , "Fall17_17Nov2017F_V32_DATA"},
     {"2018A" , "Autumn18_RunA_V19_DATA"},
     {"2018B" , "Autumn18_RunB_V19_DATA"},
     {"2018C" , "Autumn18_RunC_V19_DATA"},
     {"2018D" , "Autumn18_RunD_V19_DATA"}
   };

   static std::map<std::string, std::string> jerTagsMC = {
     {"2016" , "Summer16_25nsV1_MC"},
     {"2017" , "Fall17_V3_MC"},
     {"2018" , "Autumn18_V7_MC"}
   };

   std::string jertag = jerTagsMC.begin()->second;
   std::string jectag = jecTagsMC.begin()->second;
   std::string jectagfastsim = jecTagsFastSim.begin()->second;
   std::string archivetag = archiveTagsDATA.begin()->second;


   if(isData){
     std::string delimiter = ("_Run"+year).c_str();
     unsigned int pos = outputfilename.find(delimiter.c_str()) + delimiter.length();
     runera = (year+outputfilename.substr(pos,1)).c_str();

     jectag = jecTagsDATA[runera];
     jertag = jerTagsMC[year];
     archivetag = archiveTagsDATA[year];
   }
   else{
     runera = (year+"MC").c_str();
     jertag = jerTagsMC[year];
     jectag = jecTagsMC[year];
     jectagfastsim = jecTagsFastSim[year];
   }

   try{
     jetScaleRes = JetScaleResolution((PUSPACE+"JetResDatabase/textFiles/"+jectag+"_Uncertainty_AK4PFchs.txt").c_str(),"",(PUSPACE+"JetResDatabase/textFiles/"+jertag+"_PtResolution_AK4PFchs.txt").c_str(), (PUSPACE+"JetResDatabase/textFiles/"+jertag+"_SF_AK4PFchs.txt").c_str());
   }
   catch(edm::Exception &err){
     std::cerr << "Error in setupJetCorrections (JetScaleResolution): " << err.what() << std::endl;
     std::cerr << "\tAborting Analyzer..." << std::endl;
     std::abort();
   }
   catch(...){
     std::cerr << "Error in setupJetCorrections (JetScaleResolution): unknown Exception" << std::endl;
     std::cerr << "\tAborting Analyzer..." << std::endl;
     std::abort();
    }

   try{
     // Arguments: JetRecalibrator(const std::string path, const std::string globalTag, const std::string jetFlavor, const std::string type, bool doResidualJECs, int upToLevel=3, bool calculateSeparateCorrections=false, bool calculateTypeIMETCorr=false);
     jetRecalib = JetRecalibrator((PUSPACE+"JetResDatabase/textFiles/").c_str(), jectag, "AK4PFchs", "Total", true);
     jetRecalibL1 = JetRecalibrator((PUSPACE+"JetResDatabase/textFiles/").c_str(), jectag, "AK4PFchs", "Total", false, 1, true);
   }
   catch(std::runtime_error& err){
     std::cerr << "Error in setupJetCorrections (JetRecalibrator): " << err.what() << std::endl;
     std::cerr << "\tAborting Analyzer..." << std::endl;
     std::abort();

   }
   catch(...){
     std::cerr << "Error in setupJetCorrections (JetRecalibrator): unknown Exception" << std::endl;
     std::cerr << "\tAborting Analyzer..." << std::endl;
     std::abort();
   }

 }

// --- Function that applies the latest JECs and propagates them to MET  --- //
void Analyzer::applyJetEnergyCorrections(Particle& jet, const CUTS eGenPos, const PartStats& stats, std::string year, int syst){

  if(!jet.needSyst(syst)){
    return;
  }
  else if(jet.type != PType::Jet){
    // Return if it's FatJet or something else
    jet.setOrigReco();
    return;
  }

  std::string systname = syst_names.at(syst);

  // Include the xy-shift corrections to raw MET -- this is equivalent to apply them to the type-I MET.
  if(distats["Run"].bfind("ApplyMETxyShiftCorrections")){
    _MET->applyXYshiftCorr(year, runera, totalVertices, isData, systname, syst);
  }


  // Define the deltas to be applied to MET at the end:
  float delta_x_EEnoise_rawJets = 0.0, delta_y_EEnoise_rawJets = 0.0, delta_x_EEnoise_T1Jets = 0.0, delta_y_EEnoise_T1Jets = 0.0;
  // Define the jet energy threshold below which we consider it to be unclustered energy.
  float jetUnclEnThreshold = 15.0;

  jets_jer_sfs.clear();
  jets_jetptres.clear();
  jets_datatoMCSFs_jer.clear();
  genMatchedJets.clear();
  genMatchedJets.shrink_to_fit();
  genUnmatchedJets.clear();
  genUnmatchedJets.shrink_to_fit();

  minDeltaPhiMet_formet = 9999.9;
  maxDeltaPhiMet_formet = 0.0;
  maxjetptprojonmet_plus_formet = 0.0, maxjetptprojonmet_minus_formet = 0.0;
  index_minjmetdphi_formet = -1, index_maxjmetdphi_formet = -1;
  index_maxjetptprojonmet_plus_formet = -1, index_maxjetptprojonmet_minus_formet = -1;

  for(size_t i = 0; i < jet.size(); i++){

    // Get the reconstruced 4-vector (original vector straight from the corresponding branches), these jets have applied L1L2L3 JECs
    const TLorentzVector origJetReco = jet.RecoP4(i);

    // Step 1: get the raw jet 4-momentum from the raw factor
    float jet_RawFactor = _Jet->rawFactor[i];
    float jet_Pt = origJetReco.Pt(), jet_Mass = origJetReco.M(); // L1L2L3 corrected.
    float jet_rawPt = jet_Pt * (1.0 - jet_RawFactor), jet_rawMass = jet_Mass * (1.0 - jet_RawFactor); // raw jet pt and mass (no corrections)
    TLorentzVector rawJetP4(0,0,0,0);
    rawJetP4.SetPtEtaPhiM(jet_rawPt, origJetReco.Eta(), origJetReco.Phi(), jet_rawMass); // This is our raw jet 4-momentum

    // Step 2: check if there is any muon overlapping with this jet. If so, remove the muon(s) 4-momentum(a) from the raw jet 4-momentum
    TLorentzVector rawJetP4_noMuon = rawJetP4;
    TLorentzVector muonsP4(0,0,0,0);

    if( (_Jet->matchingMuonIdx1[i] > -1) && (_Muon->isGlobal[_Jet->matchingMuonIdx1[i]] == true) ){
      rawJetP4_noMuon = rawJetP4_noMuon - _Muon->p4(_Jet->matchingMuonIdx1[i]);
      muonsP4 += _Muon->p4(_Jet->matchingMuonIdx1[i]);
    }

    if( (_Jet->matchingMuonIdx2[i] > -1) && (_Muon->isGlobal[_Jet->matchingMuonIdx2[i]] == true) ){
      rawJetP4_noMuon = rawJetP4_noMuon - _Muon->p4(_Jet->matchingMuonIdx2[i]);
      muonsP4 += _Muon->p4(_Jet->matchingMuonIdx2[i]);
    }

    // Step 3: check if the jet pt of the corrected rawJetP4 with L1L2L3 corrections is above the unclustered energy threshold.

    // Step 3.1: get the JEC L1L2L3 and JEC L1 factors for this jet (L1 to be used later)
    double jecL1L2L3 = jetRecalib.getCorrection(rawJetP4_noMuon, _Jet->area[i], jec_rho);
    double jecL1 = jetRecalibL1.getCorrection(rawJetP4_noMuon, _Jet->area[i], jec_rho);

    // Step 3.2: check if is above the unclustered threshold. If it is, then apply the correction factor, otherwise, don't apply it.
    TLorentzVector jetL1L2L3_noMuonP4(0,0,0,0);
    TLorentzVector jetL1_noMuonP4(0,0,0,0);

    if(year.compare("2017") != 0 ){
      // For 2016 and 2018, apply the corrections as usual.
      if( jecL1L2L3 * rawJetP4_noMuon.Pt() > jetUnclEnThreshold ){
        jetL1L2L3_noMuonP4 = jetRecalib.correctedP4(rawJetP4_noMuon, jecL1L2L3); // L1L2L3 correction.
        jetL1_noMuonP4 = jetRecalibL1.correctedP4(rawJetP4_noMuon, jecL1); // L1 correction.
      } else {
        jetL1L2L3_noMuonP4 = rawJetP4_noMuon; // No correction.
        jetL1_noMuonP4 = rawJetP4_noMuon; // No correction.
      }
    } else {
      // This step is only needed for v2 MET in 2017, when different JECs are applied compared to the nanoAOD production.
      // only correct the non-mu fraction of the jet. If the corrected pt > 15 GeV (unclEnThreshold), use the corrected jet, otherwise use raw
      if(jecL1L2L3 * rawJetP4_noMuon.Pt() > jetUnclEnThreshold){
        jetL1L2L3_noMuonP4 = jetRecalib.correctedP4(rawJetP4_noMuon, jecL1L2L3); // L1L2L3 correction.
      } else {
        jetL1L2L3_noMuonP4 = rawJetP4_noMuon; // No correction.
      }

      if(jecL1 * rawJetP4_noMuon.Pt() > jetUnclEnThreshold){
        jetL1_noMuonP4 = jetRecalib.correctedP4(rawJetP4_noMuon, jecL1); // L1 correction.
      } else {
        jetL1_noMuonP4 = rawJetP4_noMuon; // No correction.
      }
    }

    // Step 4 (optional): Apply the JER corrections if desired to MC
    // Define the JER scale factors:
    double jer_sf_nom = 1.0, jer_shift = 1.0;
    // The JER corrections should be applied on top of the JECs. Then, we must use the jec corrected pt/mass of the jet without the muon!
    TLorentzVector jetL1L2L3_jerNom_noMuonP4(0,0,0,0), jetL1L2L3_jerShifted_noMuonP4(0,0,0,0);

    jetL1L2L3_jerNom_noMuonP4.SetPtEtaPhiM( jetL1L2L3_noMuonP4.Pt()*jer_sf_nom , jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jetL1L2L3_noMuonP4.M()*jer_sf_nom );
    jetL1L2L3_jerShifted_noMuonP4.SetPtEtaPhiM( jetL1L2L3_noMuonP4.Pt()*jer_shift , jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jetL1L2L3_noMuonP4.M()*jer_shift );

    // Define some helpful variables.
    float jet_pt_nomu_nom = jetL1L2L3_jerNom_noMuonP4.Pt(), jet_mass_nomu_nom = jetL1L2L3_jerNom_noMuonP4.M();
    float jet_pt_nomu_jerShifted = jetL1L2L3_jerShifted_noMuonP4.Pt(), jet_mass_nomu_jerShifted = jetL1L2L3_jerShifted_noMuonP4.M();

    // Verify that this is done (1) for MC, (2) if the smearing is turned on, (3) if the corrected jet Pt without the muon is above the unclustered energy threshold.
    if(isData && (jetL1L2L3_noMuonP4.Pt() > jetUnclEnThreshold) ){
      // study separation of jets from Met
      if( (year.compare("2016") == 0) || (year.compare("2018") == 0)){
        float deltaPhiJetMet_data = normPhi(rawJetP4_noMuon.Phi() - _MET->phi());
        float cosinedPhiJetMet_data = cos(deltaPhiJetMet_data);
        float jetptmetproj_p_data = 0.0, jetptmetproj_m_data = 0.0;

        if(cosinedPhiJetMet_data < 0.0){
          jetptmetproj_m_data = rawJetP4_noMuon.Pt() * cosinedPhiJetMet_data;
        } else if(cosinedPhiJetMet_data > 0.0){
          jetptmetproj_p_data = rawJetP4_noMuon.Pt() * cosinedPhiJetMet_data;
        }

        if(abs(deltaPhiJetMet_data) < minDeltaPhiMet_formet){
          minDeltaPhiMet_formet = abs(deltaPhiJetMet_data);
          index_minjmetdphi_formet = i;
        }

        if(abs(deltaPhiJetMet_data) > maxDeltaPhiMet_formet){
          maxDeltaPhiMet_formet = abs(deltaPhiJetMet_data);
          index_maxjmetdphi_formet = i;
        }

        if(abs(jetptmetproj_m_data) > maxjetptprojonmet_minus_formet){
          maxjetptprojonmet_minus_formet = abs(jetptmetproj_m_data);
          index_maxjetptprojonmet_minus_formet = i;
        }

        if(abs(jetptmetproj_p_data) > maxjetptprojonmet_plus_formet){
          maxjetptprojonmet_plus_formet = abs(jetptmetproj_p_data);
          index_maxjetptprojonmet_plus_formet = i;
        }
      } else if ( year.compare("2017") == 0 ){
        if( ! ((abs(rawJetP4_noMuon.Eta()) > 2.65 && abs(rawJetP4_noMuon.Eta()) < 3.14 ) && (rawJetP4_noMuon.Pt() < 50.0) ) ){
          float deltaPhiJetMet_data = normPhi(rawJetP4_noMuon.Phi() - _MET->phi());
          float cosinedPhiJetMet_data = cos(deltaPhiJetMet_data);
          float jetptmetproj_p_data = 0.0, jetptmetproj_m_data = 0.0;

          if(cosinedPhiJetMet_data < 0.0){
            jetptmetproj_m_data = rawJetP4_noMuon.Pt() * cosinedPhiJetMet_data;
          } else if(cosinedPhiJetMet_data > 0.0){
            jetptmetproj_p_data = rawJetP4_noMuon.Pt() * cosinedPhiJetMet_data;
          }

          if(abs(deltaPhiJetMet_data) < minDeltaPhiMet_formet){
            minDeltaPhiMet_formet = abs(deltaPhiJetMet_data);
            index_minjmetdphi_formet = i;
          }

          if(abs(deltaPhiJetMet_data) > maxDeltaPhiMet_formet){
            maxDeltaPhiMet_formet = abs(deltaPhiJetMet_data);
            index_maxjmetdphi_formet = i;
          }

          if(abs(jetptmetproj_m_data) > maxjetptprojonmet_minus_formet){
            maxjetptprojonmet_minus_formet = abs(jetptmetproj_m_data);
            index_maxjetptprojonmet_minus_formet = i;
          }

          if(abs(jetptmetproj_p_data) > maxjetptprojonmet_plus_formet){
            maxjetptprojonmet_plus_formet = abs(jetptmetproj_p_data);
            index_maxjetptprojonmet_plus_formet = i;
          }
        }
      }
    }

    if( (!isData) && (stats.bfind("SmearTheJet")) && (jetL1L2L3_noMuonP4.Pt() > jetUnclEnThreshold) ){

      // Check if there are any overlaps only with muons. Treat any electrons or hadronic taus which show up in the jet vector as jets.
      // Jet overlap removal in the config files takes care of them later by removing them of the jet vector.
      bool jetlepmatch = false;
      // Check that this jet doesn't match a lepton at gen-level. This will make sure that the reco jet is a true jet.
      if(JetMatchesLepton(*_Muon, origJetReco, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon)){
        // Investigate if this is a muon coming from a b-jet by checking the flavor of the gen-level jet:
        if(_Jet->genJetIdx[i] != -1){
          // Get the hadron flavor:
          int jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[i]);
          if( abs(jetHadronFlavor) != 5 || abs(_GenJet->genPartonFlavor[_Jet->genJetIdx[i]]) != 5){ // if not a b-jet, then lepton match is true.
            jetlepmatch = true;
          }
        } else {
          jetlepmatch = true;
        }
      }

      // study separation of jets from Met
      if(!jetlepmatch){

        if( (year.compare("2016") == 0) || (year.compare("2018") == 0)){
          float deltaPhiJetMet = normPhi(rawJetP4_noMuon.Phi() - _MET->phi());
          float cosinedPhiJetMet = cos(deltaPhiJetMet);
          float jetptmetproj_p = 0.0, jetptmetproj_m = 0.0;

          if(cosinedPhiJetMet < 0.0){
            jetptmetproj_m = rawJetP4_noMuon.Pt() * cosinedPhiJetMet;
          } else if(cosinedPhiJetMet > 0.0){
            jetptmetproj_p = rawJetP4_noMuon.Pt() * cosinedPhiJetMet;
          }

          if(abs(deltaPhiJetMet) < minDeltaPhiMet_formet){
            minDeltaPhiMet_formet = abs(deltaPhiJetMet);
            index_minjmetdphi_formet = i;
          }

          if(abs(deltaPhiJetMet) > maxDeltaPhiMet_formet){
            maxDeltaPhiMet_formet = abs(deltaPhiJetMet);
            index_maxjmetdphi_formet = i;
          }

          if(abs(jetptmetproj_m) > maxjetptprojonmet_minus_formet){
            maxjetptprojonmet_minus_formet = abs(jetptmetproj_m);
            index_maxjetptprojonmet_minus_formet = i;
          }

          if(abs(jetptmetproj_p) > maxjetptprojonmet_plus_formet){
            maxjetptprojonmet_plus_formet = abs(jetptmetproj_p);
            index_maxjetptprojonmet_plus_formet = i;
          }
        } else if ( year.compare("2017") == 0 ){
          if( ! ((abs(rawJetP4_noMuon.Eta()) > 2.65 && abs(rawJetP4_noMuon.Eta()) < 3.14 ) && (rawJetP4_noMuon.Pt() < 50.0) ) ){
            float deltaPhiJetMet = normPhi(rawJetP4_noMuon.Phi() - _MET->phi());
            float cosinedPhiJetMet = cos(deltaPhiJetMet);
            float jetptmetproj_p = 0.0, jetptmetproj_m = 0.0;

            if(cosinedPhiJetMet < 0.0){
              jetptmetproj_m = rawJetP4_noMuon.Pt() * cosinedPhiJetMet;
            } else if(cosinedPhiJetMet > 0.0){
              jetptmetproj_p = rawJetP4_noMuon.Pt() * cosinedPhiJetMet;
            }

            if(abs(deltaPhiJetMet) < minDeltaPhiMet_formet){
              minDeltaPhiMet_formet = abs(deltaPhiJetMet);
              index_minjmetdphi_formet = i;
            }

            if(abs(deltaPhiJetMet) > maxDeltaPhiMet_formet){
              maxDeltaPhiMet_formet = abs(deltaPhiJetMet);
              index_maxjmetdphi_formet = i;
            }

            if(abs(jetptmetproj_m) > maxjetptprojonmet_minus_formet){
              maxjetptprojonmet_minus_formet = abs(jetptmetproj_m);
              index_maxjetptprojonmet_minus_formet = i;
            }

            if(abs(jetptmetproj_p) > maxjetptprojonmet_plus_formet){
              maxjetptprojonmet_plus_formet = abs(jetptmetproj_p);
              index_maxjetptprojonmet_plus_formet = i;
            }
          }
        }
      }

      float genJetMatchDR = 0.0;
      try{
        genJetMatchDR = jet.pstats["Smear"].dmap.at("GenMatchingDeltaR");
      }
      catch(std::out_of_range& err){
        if(jet.type == PType::Jet) genJetMatchDR = 0.2;
        else if(jet.type == PType::FatJet) genJetMatchDR = 0.4;
      }

      // Find the gen-level jet that matches this reco jet (the original jet).
      TLorentzVector genJet = matchJetToGen(jetL1L2L3_noMuonP4, genJetMatchDR, eGenPos, stats.bfind("ResolutionMatching"));

      if(genJet != TLorentzVector(0,0,0,0)){
        genMatchedJets.push_back(i);
      } else {
        genUnmatchedJets.push_back(i);
      }
      // Use the jet without the muon momentum to retrieve the appropriate JER scale factors.
      // Save the three JER scale factors (nominal, down and up) as a vector in a vector: 0 - nominal, 1 - down, 2 - up
      std::vector<float> jet_jer_sf = jetScaleRes.GetSmearValsPtSF(jetL1L2L3_noMuonP4, genJet, jec_rho);
      jets_jer_sfs[i] = jet_jer_sf;

      jets_jetptres[i] = jetScaleRes.getRelativePtResMC(jetL1L2L3_noMuonP4, jec_rho);

      if(systname == "orig"){
        jets_datatoMCSFs_jer[i] = jetScaleRes.getDataToMCCoreResSF(jetL1L2L3_noMuonP4, 0);
      } else if(systname == "Jet_Res_Up"){
        jets_datatoMCSFs_jer[i] = jetScaleRes.getDataToMCCoreResSF(jetL1L2L3_noMuonP4, 2);
      } else if(systname == "Jet_Res_Down"){
        jets_datatoMCSFs_jer[i] = jetScaleRes.getDataToMCCoreResSF(jetL1L2L3_noMuonP4, 1);
      }

      bool genjetmatch = false, passtightPUjetID = false, smearunmatchedjet = true, forwardjetsmearing = true;

      if( stats.bfind("ModifiedPUsmearing") ){

        genjetmatch = genJet != TLorentzVector(0,0,0,0);
        passtightPUjetID = _Jet->getPileupJetID(i, 0); // i - jet index, 0 = bit0 (tight ID). If true, it passes tight PU jet ID, otherwise, it fails PU jet ID.
        if (jetL1L2L3_noMuonP4.Pt() <= 50.0) smearunmatchedjet = !genjetmatch && passtightPUjetID ? true : false;

        if( (genjetmatch || smearunmatchedjet) && !jetlepmatch){

          // Update the resolution scale factors as well as the 4 vectors for the jet momentum.
          jer_sf_nom = jet_jer_sf.at(0); // 0 is the nominal value

          // Correct the nominal mass and pt for this jet.
          jet_pt_nomu_nom = jetL1L2L3_noMuonP4.Pt() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.Pt() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.Pt() * jer_sf_nom;
          jet_mass_nomu_nom = jetL1L2L3_noMuonP4.M() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.M() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.M() * jer_sf_nom;

          if(systname == "orig"){ // This corresponds to the nominal values: set the scale factor; if smearing, update jet_pt_nom and jet_mass_nom
            jer_shift = jer_sf_nom;
            jet_pt_nomu_jerShifted = jet_pt_nomu_nom;
            jet_mass_nomu_jerShifted = jet_mass_nomu_nom;

          }else if(systname.find("_Res_") != std::string::npos){
            if(systname == "Jet_Res_Up"){
              jer_shift = jet_jer_sf.at(2); // i is the jet index, 2 is the up value
            }else if(systname == "Jet_Res_Down"){
              jer_shift = jet_jer_sf.at(1); // i is the jet index, 1 is the down value
            }
            jet_pt_nomu_jerShifted = jetL1L2L3_noMuonP4.Pt() * jer_shift;
            jet_mass_nomu_jerShifted = jetL1L2L3_noMuonP4.M() * jer_shift;
          }
        }

      } else if ( stats.bfind("ModifiedForwardSmearing") ){

        if( (jetL1L2L3_noMuonP4.Pt() <= 50.0) && ( abs(jetL1L2L3_noMuonP4.Eta()) > 2.5 ) ) forwardjetsmearing = false;

        if( forwardjetsmearing && !jetlepmatch){

          // Update the resolution scale factors as well as the 4 vectors for the jet momentum.
          jer_sf_nom = jet_jer_sf.at(0); // 0 is the nominal value

          // Correct the nominal mass and pt for this jet.
          jet_pt_nomu_nom = jetL1L2L3_noMuonP4.Pt() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.Pt() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.Pt() * jer_sf_nom;
          jet_mass_nomu_nom = jetL1L2L3_noMuonP4.M() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.M() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.M() * jer_sf_nom;

          if(systname == "orig"){ // This corresponds to the nominal values: set the scale factor; if smearing, update jet_pt_nom and jet_mass_nom
            jer_shift = jer_sf_nom;
            jet_pt_nomu_jerShifted = jet_pt_nomu_nom;
            jet_mass_nomu_jerShifted = jet_mass_nomu_nom;

          }else if(systname.find("_Res_") != std::string::npos){
            if(systname == "Jet_Res_Up"){
              jer_shift = jet_jer_sf.at(2); // i is the jet index, 2 is the up value
            }else if(systname == "Jet_Res_Down"){
              jer_shift = jet_jer_sf.at(1); // i is the jet index, 1 is the down value
            }
            jet_pt_nomu_jerShifted = jetL1L2L3_noMuonP4.Pt() * jer_shift;
            jet_mass_nomu_jerShifted = jetL1L2L3_noMuonP4.M() * jer_shift;
          }
        }

      } else if ( stats.bfind("CombinedModifiedSmearing") ){

        if( abs( jetL1L2L3_noMuonP4.Eta() ) <= 2.5 ){ // central jets

          genjetmatch = genJet != TLorentzVector(0,0,0,0);
          passtightPUjetID = _Jet->getPileupJetID(i, 0); // i - jet index, 0 = bit0 (tight ID). If true, it passes tight PU jet ID, otherwise, it fails PU jet ID.
          if (jetL1L2L3_noMuonP4.Pt() <= 50.0) smearunmatchedjet = !genjetmatch && passtightPUjetID ? true : false;

          if( (genjetmatch || smearunmatchedjet) && !jetlepmatch){

            // Update the resolution scale factors as well as the 4 vectors for the jet momentum.
            jer_sf_nom = jet_jer_sf.at(0); // 0 is the nominal value

            // Correct the nominal mass and pt for this jet.
            jet_pt_nomu_nom = jetL1L2L3_noMuonP4.Pt() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.Pt() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.Pt() * jer_sf_nom;
            jet_mass_nomu_nom = jetL1L2L3_noMuonP4.M() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.M() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.M() * jer_sf_nom;

            if(systname == "orig"){ // This corresponds to the nominal values: set the scale factor; if smearing, update jet_pt_nom and jet_mass_nom
              jer_shift = jer_sf_nom;
              jet_pt_nomu_jerShifted = jet_pt_nomu_nom;
              jet_mass_nomu_jerShifted = jet_mass_nomu_nom;

            }else if(systname.find("_Res_") != std::string::npos){
              if(systname == "Jet_Res_Up"){
                jer_shift = jet_jer_sf.at(2); // i is the jet index, 2 is the up value
              }else if(systname == "Jet_Res_Down"){
                jer_shift = jet_jer_sf.at(1); // i is the jet index, 1 is the down value
              }
              jet_pt_nomu_jerShifted = jetL1L2L3_noMuonP4.Pt() * jer_shift;
              jet_mass_nomu_jerShifted = jetL1L2L3_noMuonP4.M() * jer_shift;
            }
          }

        } else if( abs( jetL1L2L3_noMuonP4.Eta() ) > 2.5 ) { // forward jets

          if( jetL1L2L3_noMuonP4.Pt() <= 50.0 ) forwardjetsmearing = false;

          if( forwardjetsmearing && !jetlepmatch){

            // Update the resolution scale factors as well as the 4 vectors for the jet momentum.
            jer_sf_nom = jet_jer_sf.at(0); // 0 is the nominal value

            // Correct the nominal mass and pt for this jet.
            jet_pt_nomu_nom = jetL1L2L3_noMuonP4.Pt() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.Pt() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.Pt() * jer_sf_nom;
            jet_mass_nomu_nom = jetL1L2L3_noMuonP4.M() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.M() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.M() * jer_sf_nom;

            if(systname == "orig"){ // This corresponds to the nominal values: set the scale factor; if smearing, update jet_pt_nom and jet_mass_nom
              jer_shift = jer_sf_nom;
              jet_pt_nomu_jerShifted = jet_pt_nomu_nom;
              jet_mass_nomu_jerShifted = jet_mass_nomu_nom;

            }else if(systname.find("_Res_") != std::string::npos){
              if(systname == "Jet_Res_Up"){
                jer_shift = jet_jer_sf.at(2); // i is the jet index, 2 is the up value
              }else if(systname == "Jet_Res_Down"){
                jer_shift = jet_jer_sf.at(1); // i is the jet index, 1 is the down value
              }
              jet_pt_nomu_jerShifted = jetL1L2L3_noMuonP4.Pt() * jer_shift;
              jet_mass_nomu_jerShifted = jetL1L2L3_noMuonP4.M() * jer_shift;
            }
          }

        }

      } else {

        // Update the resolution scale factors as well as the 4 vectors for the jet momentum.
        jer_sf_nom = jet_jer_sf.at(0); // 0 is the nominal value

        jet_pt_nomu_nom = jetL1L2L3_noMuonP4.Pt() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.Pt() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.Pt() * jer_sf_nom;
        jet_mass_nomu_nom = jetL1L2L3_noMuonP4.M() * jer_sf_nom > 0.0 ? jetL1L2L3_noMuonP4.M() * jer_sf_nom : -1.0 * jetL1L2L3_noMuonP4.M() * jer_sf_nom;

        if(!jetlepmatch){ // Do the smearing only if there is no overlap with muons, this block will only change the sfs and jet pt/mass.

          if(systname == "orig"){ // This corresponds to the nominal values; set the scale factor; if smearing, update jet_pt_nom and jet_mass_nom
            jer_shift = jer_sf_nom;
            jet_pt_nomu_jerShifted = jet_pt_nomu_nom;
            jet_mass_nomu_jerShifted = jet_mass_nomu_nom;

          } else if(systname.find("_Res_") != std::string::npos){
            if(systname == "Jet_Res_Up"){
              jer_shift = jet_jer_sf.at(2); // i is the jet index, 2 is the up value
            }else if(systname == "Jet_Res_Down"){
              jer_shift = jet_jer_sf.at(1); // i is the jet index, 1 is the down value
            }
            jet_pt_nomu_jerShifted = jetL1L2L3_noMuonP4.Pt() * jer_shift;
            jet_mass_nomu_jerShifted = jetL1L2L3_noMuonP4.M() * jer_shift;
          }

        }

      }

    }

    //----------
    // Update the pt and mass values for the nominal jet P4:
    jetL1L2L3_jerNom_noMuonP4.SetPtEtaPhiM( jet_pt_nomu_nom , jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jet_mass_nomu_nom );
    // Update the 4-vector for the jet with shifted resolution.
    jetL1L2L3_jerShifted_noMuonP4.SetPtEtaPhiM( jet_pt_nomu_jerShifted , jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jet_mass_nomu_jerShifted );

    // After correcting and smearing the clustered energy, we can now add back in the muon momenta, if any.
    TLorentzVector jetL1L2L3_jerNomP4 = jetL1L2L3_jerNom_noMuonP4 + muonsP4;
    TLorentzVector jetL1L2L3_jerShiftedP4 = jetL1L2L3_jerShifted_noMuonP4 + muonsP4;
    TLorentzVector jetL1_P4 = jetL1_noMuonP4 + muonsP4;

    // Evaluate JES uncertainties. This must be done on top of the nominal JECs and JERs (if any), not considering the muon momenta.
    TLorentzVector jetL1L2L3_jer_noMuon_jesShiftedP4(0,0,0,0), jetL1L2L3_T1_noMuon_jesShiftedP4(0,0,0,0);
    double jes_delta = 0.0, jes_delta_t1 = 0.0;

    float jet_pt_jesShifted = jetL1L2L3_jerNom_noMuonP4.Pt() * (1.0 + jes_delta), jet_mass_jesShifted = jetL1L2L3_jerNom_noMuonP4.M() * (1.0 + jes_delta);
    float jet_pt_jesShiftedT1 = jetL1L2L3_noMuonP4.Pt() * (1.0 + jes_delta_t1), jet_mass_jesShiftedT1 = jetL1L2L3_noMuonP4.M() * (1.0 + jes_delta_t1);

    jetL1L2L3_jer_noMuon_jesShiftedP4.SetPtEtaPhiM( jetL1L2L3_jerNom_noMuonP4.Pt() * (1.0 + jes_delta), jetL1L2L3_jerNom_noMuonP4.Eta(), jetL1L2L3_jerNom_noMuonP4.Phi(), jetL1L2L3_jerNom_noMuonP4.M() * (1.0 + jes_delta) );
    jetL1L2L3_T1_noMuon_jesShiftedP4.SetPtEtaPhiM( jetL1L2L3_noMuonP4.Pt() * (1.0 + jes_delta_t1), jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jetL1L2L3_noMuonP4.M() * (1.0 + jes_delta_t1) );

    if(systname.find("_Scale_") != std::string::npos){

      // Here we will be using the mass and pt that were obtained after applying JER corrections
      jes_delta = jetScaleRes.GetScaleDelta(jetL1L2L3_jerNom_noMuonP4.Pt(), jetL1L2L3_jerNom_noMuonP4.Eta());
      jes_delta_t1 = jetScaleRes.GetScaleDelta(jetL1L2L3_noMuonP4.Pt(), jetL1L2L3_noMuonP4.Eta());

      // JES applied for systematics both in data and MC.
      if(systname == "Jet_Scale_Up"){
        jet_pt_jesShifted = jetL1L2L3_jerNom_noMuonP4.Pt() * (1.0 + jes_delta);
        jet_mass_jesShifted = jetL1L2L3_jerNom_noMuonP4.M() * (1.0 + jes_delta);

        // If no smearing is applied, just re-do JES variations for T1 MET
        jet_pt_jesShiftedT1 = jetL1L2L3_noMuonP4.Pt() * (1.0 + jes_delta_t1);
        jet_mass_jesShiftedT1 = jetL1L2L3_noMuonP4.M() * (1.0 + jes_delta_t1);

      }
      else if(systname == "Jet_Scale_Down"){
        jet_pt_jesShifted = jetL1L2L3_jerNom_noMuonP4.Pt() * (1.0 - jes_delta);
        jet_mass_jesShifted = jetL1L2L3_jerNom_noMuonP4.M() * (1.0 - jes_delta);

        // If no smearing is applied, just re-do JES variations for T1 MET
        jet_pt_jesShiftedT1 = jetL1L2L3_noMuonP4.Pt() * (1.0 - jes_delta_t1);
        jet_mass_jesShiftedT1 = jetL1L2L3_noMuonP4.M() * (1.0 - jes_delta_t1);
      }
    }

    // Update the shifted JES vectors accordingly
    jetL1L2L3_jer_noMuon_jesShiftedP4.SetPtEtaPhiM( jet_pt_jesShifted, jetL1L2L3_jerNom_noMuonP4.Eta(), jetL1L2L3_jerNom_noMuonP4.Phi(), jet_mass_jesShifted );
    jetL1L2L3_T1_noMuon_jesShiftedP4.SetPtEtaPhiM( jet_pt_jesShiftedT1, jetL1L2L3_noMuonP4.Eta(), jetL1L2L3_noMuonP4.Phi(), jet_mass_jesShiftedT1 );

    // Add back in the muon momenta, if any.
    TLorentzVector jetL1L2L3_jer_jesShiftedP4 = jetL1L2L3_jer_noMuon_jesShiftedP4 + muonsP4;
    TLorentzVector jetL1L2L3_T1_jesShiftedP4 = jetL1L2L3_T1_noMuon_jesShiftedP4 + muonsP4;

    // Shift the jet 4-momentum accordingly:
    if(isData){
    	// Here, jet_pt_jerShifted and jet_mass_jerShifted are unchanged w.r.t. the original values. We need to do this in order to keep these jets in the _Jet vector.
    	systematics.shiftJet(jet, jetL1L2L3_jerShiftedP4, systname, syst);
    }
    else if(!isData){
    	if(stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos){
    		// Correct the jet 4-momentum according to the systematic applied for JES.
      	systematics.shiftJet(jet, jetL1L2L3_jer_jesShiftedP4, systname, syst);
    	} else if(!stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos){
        // Correct the jet 4-momentum according to the systematic applied for JES.
        systematics.shiftJet(jet, jetL1L2L3_T1_jesShiftedP4, systname, syst);
      } else {
    		// Correct the jet 4-momentum according to the systematic applied for JER
        systematics.shiftJet(jet, jetL1L2L3_jerShiftedP4, systname, syst);
    	}
    }

    // --------------------- Propagation of JER/JES to MET --------------------- //

    double jetTotalEmEF = _Jet->neutralEmEnergyFraction[i] + _Jet->chargedEmEnergyFraction[i];

    // Propagate JER and JES corrections and uncertainties to MET.
    // Only propagate JECs to MET if the corrected pt without the muon is above the unclustered energy threshold.

    if(jetL1L2L3_noMuonP4.Pt() > jetUnclEnThreshold && jetTotalEmEF < 0.9){

      if( ( ( year.compare("2016") == 0 ) || (year.compare("2018") == 0) ) ){
        // Do the normal propagation of JEC & JER to 2016 and 2018 and those jets outside the problematic EE noise jet eta and pt regions for 2017.
        if(systname.find("orig") != std::string::npos ){
          _MET->propagateJetEnergyCorr(jetL1L2L3_jerNom_noMuonP4, jetL1L2L3_jerNom_noMuonP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
        } else if( !stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos ){
          _MET->propagateJetEnergyCorr(jetL1L2L3_T1_noMuon_jesShiftedP4, jetL1L2L3_T1_noMuon_jesShiftedP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
        } else if( stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos ){
          _MET->propagateJetEnergyCorr(jetL1L2L3_jer_noMuon_jesShiftedP4, jetL1L2L3_jer_noMuon_jesShiftedP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
        } else if( stats.bfind("SmearTheJet") && systname.find("_Res_") != std::string::npos ){
          _MET->propagateJetEnergyCorr(jetL1L2L3_jerShifted_noMuonP4, jetL1L2L3_jerShifted_noMuonP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
        }

      } else if( year.compare("2017") == 0 ){
        if( ! ((abs(rawJetP4_noMuon.Eta()) > 2.65 && abs(rawJetP4_noMuon.Eta()) < 3.14 ) && (rawJetP4_noMuon.Pt() < 50.0) ) ){

          // Do the normal propagation of JEC & JER to 2016 and 2018 and those jets outside the problematic EE noise jet eta and pt regions for 2017.
          if(systname.find("orig") != std::string::npos ){
            _MET->propagateJetEnergyCorr(jetL1L2L3_jerNom_noMuonP4, jetL1L2L3_jerNom_noMuonP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
          } else if( !stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos ){
            _MET->propagateJetEnergyCorr(jetL1L2L3_T1_noMuon_jesShiftedP4, jetL1L2L3_T1_noMuon_jesShiftedP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
          } else if( stats.bfind("SmearTheJet") && systname.find("_Scale_") != std::string::npos ){
            _MET->propagateJetEnergyCorr(jetL1L2L3_jer_noMuon_jesShiftedP4, jetL1L2L3_jer_noMuon_jesShiftedP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
          } else if( stats.bfind("SmearTheJet") && systname.find("_Res_") != std::string::npos ){
            _MET->propagateJetEnergyCorr(jetL1L2L3_jerShifted_noMuonP4, jetL1L2L3_jerShifted_noMuonP4.Pt(), jetL1_noMuonP4.Pt(), systname, syst);
          }
        } else if( (abs(rawJetP4_noMuon.Eta()) > 2.65 && abs(rawJetP4_noMuon.Eta()) < 3.14 ) && rawJetP4_noMuon.Pt() < 50.0 ){

          // get the delta for removing raw jets in the EE region from the raw MET
          delta_x_EEnoise_rawJets += rawJetP4_noMuon.Pt() * cos(rawJetP4_noMuon.Phi());
          delta_y_EEnoise_rawJets += rawJetP4_noMuon.Pt() * sin(rawJetP4_noMuon.Phi());

          // Get the delta for removing L1L2L3-L1 corrected jets in the EE region from the default MET branch
          // Take into account if the jets are smeared in resolution, multiplying by jer_sf_nom
          //delta_x_EEnoise_T1Jets += (jetL1L2L3_jerNom_noMuonP4.Pt() - jetL1_noMuonP4.Pt()) * cos(jetL1L2L3_jerNom_noMuonP4.Phi()) + rawJetP4_noMuon.Pt() * cos(rawJetP4_noMuon.Phi());
          //delta_y_EEnoise_T1Jets += (jetL1L2L3_jerNom_noMuonP4.Pt() - jetL1_noMuonP4.Pt()) * sin(jetL1L2L3_jerNom_noMuonP4.Phi()) + rawJetP4_noMuon.Pt() * sin(rawJetP4_noMuon.Phi());
          delta_x_EEnoise_T1Jets += (jetL1L2L3_noMuonP4.Pt() - jetL1_noMuonP4.Pt()) * cos(jetL1L2L3_noMuonP4.Phi()) + rawJetP4_noMuon.Pt() * cos(rawJetP4_noMuon.Phi());
          delta_y_EEnoise_T1Jets += (jetL1L2L3_noMuonP4.Pt() - jetL1_noMuonP4.Pt()) * sin(jetL1L2L3_noMuonP4.Phi()) + rawJetP4_noMuon.Pt() * sin(rawJetP4_noMuon.Phi());

        }
      }
    }

  }

  // Propagate "unclustered energy" uncertainty to MET to finalize the 2017 MET recipe v2
  if(year.compare("2017") == 0){
    // Remove the L1L2L3-L1 corrected jets in the EE region from the default MET branch
    _MET->removeEEnoiseUnclEnergy(delta_x_EEnoise_T1Jets, delta_y_EEnoise_T1Jets, delta_x_EEnoise_rawJets, delta_y_EEnoise_rawJets, systname, syst);
  }

}

// --- Old function that smeares the jet energy resolution: this is done only in MC to improve the agreement between data and MC --- //

void Analyzer::smearJetRes(Particle& jet, const CUTS eGenPos, const PartStats& stats, int syst) {
  //at the moment
  if(isData || jet.type != PType::Jet ){
    // If it is data or not a reco jet collection, then return the original values for the 4-momenta of all particles.
     jet.setOrigReco();
     return;
  }
  else if(!jet.needSyst(syst)){
     // If this function is called but not needed, just return.
     return;
   }


  std::string systname = syst_names.at(syst);

  double genJetMatchDR = 0.0;

   try{
     genJetMatchDR = jet.pstats["Smear"].dmap.at("GenMatchingDeltaR");
   }
   catch(std::out_of_range& err){
       if(jet.type == PType::Jet) genJetMatchDR = 0.4;
       else if(jet.type == PType::FatJet) genJetMatchDR = 0.8;
   }

  for(size_t i=0; i< jet.size(); i++) {

    TLorentzVector jetReco = jet.RecoP4(i);

     // Check that this jet doesn't match a lepton at gen-level. This will make sure that the reco jet is a true jet.

    if(JetMatchesLepton(*_Muon, jetReco, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetReco, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGHadTau) ||
       JetMatchesLepton(*_Electron, jetReco,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){
      jet.addP4Syst(jetReco,syst);
      continue;
    }

    // Define the JER scale factors:
     double jer_shift = 1.0;

     // Define the new jet pt and mass variables to be updated after applying the JER corrections.
     double jet_pt_jerShifted = jetReco.Pt() * jer_shift, jet_mass_jerShifted = jetReco.M() * jer_shift;

     // Find the gen-level jet that matches this reco jet.
    TLorentzVector genJet = matchJetToGen(jetReco, genJetMatchDR, eGenPos, false);

     if(systname == "orig" && stats.bfind("SmearTheJet")){ // This corresponds to the nominal values

       // Get the scale factors:
    	jer_shift = jetScaleRes.GetSmearValsPtSF(jetReco, genJet, jec_rho).at(0); // In this case, the jer_shift is the nominal jer sf.

       // If smearing, update the jet_pt_nom and jet_mass_nom:
       jet_pt_jerShifted = jetReco.Pt() * jer_shift > 0.0 ? jetReco.Pt() * jer_shift : -1.0 * jetReco.Pt() * jer_shift;
       jet_mass_jerShifted = jetReco.M() * jer_shift > 0.0 ? jetReco.M() * jer_shift : -1.0 * jetReco.M() * jer_shift;

     }else if(systname.find("_Res_") != std::string::npos){

       if(systname == "Jet_Res_Up"){
       	jer_shift = jetScaleRes.GetSmearValsPtSF(jetReco, genJet, jec_rho).at(2);
       }else if(systname == "Jet_Res_Down"){
        jer_shift = jetScaleRes.GetSmearValsPtSF(jetReco, genJet, jec_rho).at(1);
       }

       jet_pt_jerShifted = jer_shift * jetReco.Pt();
       jet_mass_jerShifted = jer_shift * jetReco.M();
     }

     // Correct the jet 4-momentum according to the systematic applied for JER
     systematics.shiftParticle(jet, jetReco, jet_pt_jerShifted, jet_mass_jerShifted, systname, syst);

   }

 }


/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {

  for(size_t j = 0; j < lepton.size(); j++) {
    if(jetV.DeltaR(lepton.RecoP4(j)) < partDeltaR && matchLeptonToGen(lepton.RecoP4(j), lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)){
      return true;
    }
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle

TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& recoLepton4Vector, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(recoLepton4Vector, stats.dmap.at("GenMatchingDeltaR"));
  }
  if(ePos == CUTS::eGHadTau){
    return matchHadTauToGen(recoLepton4Vector, stats.dmap.at("GenMatchingDeltaR"));
  }

  float minDeltaR = 999.9, tmpDeltaR = 999.9;
  int pidx = -1;

  if( !stats.bfind("UseMotherID") ){
    for(auto it : *active_part->at(ePos)) {

      tmpDeltaR = recoLepton4Vector.DeltaR(_Gen->p4(it));

      if(tmpDeltaR < minDeltaR){
        minDeltaR = tmpDeltaR;
        pidx = it;
      }
    }
  } else {
    for(auto it : *active_part->at(ePos)) {

      if( abs(_Gen->pdg_id[_Gen->genPartIdxMother[it]]) != stats.dmap.at("MotherID") ) continue;
      tmpDeltaR = recoLepton4Vector.DeltaR(_Gen->p4(it));

      if(tmpDeltaR < minDeltaR){
        minDeltaR = tmpDeltaR;
        pidx = it;
      }
    }
  }


  if(minDeltaR <= stats.dmap.at("GenMatchingDeltaR")){
    return _Gen->p4(pidx);
  }

  return TLorentzVector(0,0,0,0);
}

/*
TLorentzVector Analyzer::matchLeptonToGen(int recoLeptonIndex, const PartStats& stats, CUTS ePos) {

  if(ePos == CUTS::eGTau) {
    return matchTauToGen(recoLeptonIndex, stats.dmap.at("GenMatchingFlavor"));
  }
  if(ePos == CUTS::eGHadTau){
    return matchHadTauToGen(recoLeptonIndex, stats.dmap.at("GenMatchingDeltaR"));
  }
  for(auto it : *active_part->at(ePos)) {
    if(recoLepton4Vector.DeltaR(_Gen->p4(it)) <= stats.dmap.at("GenMatchingDeltaR")) {
      if(stats.bfind("UseMotherID") && abs(_Gen->pdg_id[_Gen->genPartIdxMother[it]]) != stats.dmap.at("MotherID")) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}
*/
///Tau specific matching function.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically

TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it !=active_part->at(CUTS::eGTau)->end();it++, i++) {
    int nu = active_part->at(CUTS::eGNuTau)->at(i);
    if(nu == -1) continue;
    // std::bitset<14> statusflags(_Gen->statusFlags[*it]);
    // std::cout << "GenPart idx = " << *it << ", id = " << _Gen->pdg_id[*it] << ", status = " << _Gen->status[*it] << ", status flags = " << statusflags << std::endl;

    genVec = _Gen->p4(*it) - _Gen->p4(nu);
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return genVec;
}

/*
TLorentzVector Analyzer::matchTauToGen(int recoTauIndex, int genPartFlavor){

  // Look for the index of the matching particle in the Gen particle collection (tau status-2, meaning not decayed)
  int matchingTauIdx = _Tau->genPartIdx[recoTauIndex];
  int tauGenPartFlavor = static_cast<unsigned>(_Tau->genPartFlav[recoTauIndex]);

  if(matchingTauIdx != -1){

    if((genPartFlav == -1) && (tauGenPartFlavor == 3 || tauGenPartFlavor == 4 || tauGenPartFlavor == 5) ){ // include all possible tau decays (tau->e decay, tau->mu decay, hadronic tau decay)
      return _Gen->p4(matchingTauIdx);
    } else if ( (genPartFlav == 0) && (tauGenPartFlavor == 0) ){
      return _Gen->p4(matchingTauIdx);
    } else {
      if(tauGenPartFlavor == genPartFlav){
        return _Gen->p4(matchingTauIdx);
      }
    }
  }

  return TLorentzVector(0,0,0,0);
}
*/
///Tau specific matching function.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchHadTauToGen(const TLorentzVector& recoTau4Vector, double recogenDeltaR) {

  for(vec_iter genhadtau_it = active_part->at(CUTS::eGHadTau)->begin(); genhadtau_it != active_part->at(CUTS::eGHadTau)->end(); genhadtau_it++){ // (genhadtau_it) is the index of the gen-level hadronic tau in the gen-hadtau vector.

    // Compare the separation between the reco and gen hadronic tau candidates. If it's greather than the requirement, continue with the next gen-tau_h candidate.
    if(recoTau4Vector.DeltaR(_GenHadTau->p4(*genhadtau_it)) > recogenDeltaR) continue;
    // If the requirement DeltaR <= recogenDeltaR, the code will get to this point. Then we want to fill a vector that only stores the matched gen taus

    if(std::find(active_part->at(CUTS::eGMatchedHadTau)->begin(), active_part->at(CUTS::eGMatchedHadTau)->end(), *genhadtau_it) == active_part->at(CUTS::eGMatchedHadTau)->end()){
    	// The if condition above will make sure to not double count the gen-taus that have already been matched. It will only fill the CUTS::eGMatchedHadTau if it's not already in there.
         active_part->at(CUTS::eGMatchedHadTau)->push_back(*genhadtau_it);
    }
    // std::cout << "Index of the mother particle in GenPart for the matched vis tau = " << _GenHadTau->genPartIdxMother[*genhadtau_it] << ", mother pdg id = " << _Gen->pdg_id[_GenHadTau->genPartIdxMother[*genhadtau_it]] << std::endl;
    active_part->at(CUTS::eGMatchedHadTau)->push_back(*genhadtau_it);

    // And we also return the gen-tau p4 vector, which we are interested in.
    return _GenHadTau->p4(*genhadtau_it);
  }
  //return genTau4Vector;
  return TLorentzVector(0,0,0,0);
}

/*
TLorentzVector Analyzer::matchHadTauToGen(const TLorentzVector& recoTau4Vector, double recogenDeltaR) {

  for(vec_iter genhadtau_it = active_part->at(CUTS::eGHadTau)->begin(); genhadtau_it != active_part->at(CUTS::eGHadTau)->end(); genhadtau_it++){ // (genhadtau_it) is the index of the gen-level hadronic tau in the gen-hadtau vector.

    // Compare the separation between the reco and gen hadronic tau candidates. If it's greather than the requirement, continue with the next gen-tau_h candidate.
    if(recoTau4Vector.DeltaR(_GenHadTau->p4(*genhadtau_it)) > recogenDeltaR) continue;
    // If the requirement DeltaR <= recogenDeltaR, the code will get to this point. Then we want to fill a vector that only stores the matched gen taus

    if(std::find(active_part->at(CUTS::eGMatchedHadTau)->begin(), active_part->at(CUTS::eGMatchedHadTau)->end(), *genhadtau_it) == active_part->at(CUTS::eGMatchedHadTau)->end()){
      // The if condition above will make sure to not double count the gen-taus that have already been matched. It will only fill the CUTS::eGMatchedHadTau if it's not already in there.
         active_part->at(CUTS::eGMatchedHadTau)->push_back(*genhadtau_it);
    }

    active_part->at(CUTS::eGMatchedHadTau)->push_back(*genhadtau_it);

    // And we also return the gen-tau p4 vector, which we are interested in.
    return _GenHadTau->p4(*genhadtau_it);
  }
  //return genTau4Vector;
  return TLorentzVector(0,0,0,0);
}
*/


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchJetToGen(const TLorentzVector& recoJet4Vector, const double& matchDeltaR, CUTS ePos, bool resmatching){

  int genjetidx = -1;
  float drMin = 999.;

  for(auto it : *active_part->at(ePos)) {
    if(resmatching){
      // Check first that the current gen jet is a candidate for matching based on resolution:
      bool is_candidate_res_matching = jetScaleRes.resolution_matching(recoJet4Vector, _GenJet->p4(it), jec_rho);
      if(!is_candidate_res_matching) continue;
    }

    // Check separation in deltaR
    float dr = recoJet4Vector.DeltaR(_GenJet->p4(it));
    if(dr < drMin){
      genjetidx = it;
      drMin = dr;
    }
  }

  if(drMin < matchDeltaR){
    return _GenJet->p4(genjetidx);
  }
  else{
    return TLorentzVector(0,0,0,0);
  }

  return TLorentzVector(0,0,0,0);

}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
int Analyzer::matchToGenPdg(const TLorentzVector& lvec, double minDR) {
  double _minDR=minDR;
  int found=-1;
  for(size_t i=0; i< _Gen->size(); i++) {

    if(lvec.DeltaR(_Gen->p4(i)) <=_minDR) {
      //only hard interaction
      if( _Gen->status[i]<10){
        found=i;
        _minDR=lvec.DeltaR(_Gen->p4(i));
      }
    }
  }
  if (found>=0){
    return _Gen->pdg_id[found];
  }
  return 0;
}


////Calculates the number of gen particles.  Based on id number and status of each particle
// The file that corresponds to PartStats& stats is Gen_info.in
void Analyzer::getGoodGen(const PartStats& stats) {

  if(! neededCuts.isPresent(CUTS::eGen)) return;

  int particle_id = 0;

  for(size_t j = 0; j < _Gen->size(); j++) {

    particle_id = abs(_Gen->pdg_id[j]);

    if( (particle_id < 5 || particle_id == 9 || particle_id == 21) && genMaper.find(particle_id) != genMaper.end() && _Gen->status[j] == genMaper.at(5)->status){
      active_part->at(genMaper.at(5)->ePos)->push_back(j);
    }
    else if(genMaper.find(particle_id) != genMaper.end() && _Gen->status[j] == genMaper.at(particle_id)->status) {

      if(particle_id == 15){

      	if(stats.bfind("DiscrTauByPtAndEta") &&  (_Gen->pt(j) < stats.pmap.at("TauPtCut").first || _Gen->pt(j) > stats.pmap.at("TauPtCut").second || abs(_Gen->eta(j)) > stats.dmap.at("TauEtaCut"))) continue;

	if(stats.bfind("DiscrTauByMotherID")){
	   int motherpart_idx = _Gen->genPartIdxMother[j];
           int mother_pid = abs(_Gen->pdg_id[motherpart_idx]);   

           if(mother_pid == particle_id){
               // std::cout << "Lepton with same ID for mother particle" << std::endl;
               int motherpart_idx_tmp = motherpart_idx;
               int mother_pid_tmp = mother_pid;            

               while(mother_pid_tmp == particle_id){
                  motherpart_idx = _Gen->genPartIdxMother[motherpart_idx_tmp];
                  mother_pid_tmp = abs(_Gen->pdg_id[motherpart_idx]);
                  motherpart_idx_tmp = motherpart_idx;
               }
             
                mother_pid = mother_pid_tmp;
                //std::cout << "Final mother ID = " << mother_pid << std::endl;
           }
           
           if( (mother_pid != stats.pmap.at("TauMotherIDs").first) && (mother_pid != stats.pmap.at("TauMotherIDs").second) ) continue;
        }
      } else if(particle_id == 11){
      	if(stats.bfind("DiscrElecByPtAndEta") &&  (_Gen->pt(j) < stats.pmap.at("ElecPtCut").first || _Gen->pt(j) > stats.pmap.at("ElecPtCut").second || abs(_Gen->eta(j)) > stats.dmap.at("ElecEtaCut"))) continue;

        if(stats.bfind("DiscrElecByMotherID")){
            int motherpart_idx = _Gen->genPartIdxMother[j];
            int mother_pid = abs(_Gen->pdg_id[motherpart_idx]);
            
            if(mother_pid == particle_id){
                //std::cout << "Lepton with same ID for mother particle" << std::endl;
                int motherpart_idx_tmp = motherpart_idx;
                int mother_pid_tmp = mother_pid;
                
                while(mother_pid_tmp == particle_id){
                   motherpart_idx = _Gen->genPartIdxMother[motherpart_idx_tmp];
                   mother_pid_tmp = abs(_Gen->pdg_id[motherpart_idx]);
                   motherpart_idx_tmp = motherpart_idx;
                } 
              
                 mother_pid = mother_pid_tmp;
                 //std::cout << "Final mother ID = " << mother_pid << std::endl;
            }
            
            if( (mother_pid != stats.pmap.at("ElecMotherIDs").first) && (mother_pid != stats.pmap.at("ElecMotherIDs").second) ) continue;
        }
      } else if(particle_id == 13){

         if(stats.bfind("DiscrMuonByPtAndEta") &&  (_Gen->pt(j) < stats.pmap.at("MuonPtCut").first || _Gen->pt(j) > stats.pmap.at("MuonPtCut").second || abs(_Gen->eta(j)) > stats.dmap.at("MuonEtaCut"))) continue;

         if(stats.bfind("DiscrTauByMotherID")){
            int motherpart_idx = _Gen->genPartIdxMother[j];
            int mother_pid = abs(_Gen->pdg_id[motherpart_idx]);
            
            if(mother_pid == particle_id){
                // std::cout << "Lepton with same ID for mother particle" << std::endl;
                int motherpart_idx_tmp = motherpart_idx;
                int mother_pid_tmp = mother_pid;
                
                while(mother_pid_tmp == particle_id){
                   motherpart_idx = _Gen->genPartIdxMother[motherpart_idx_tmp];
                   mother_pid_tmp = abs(_Gen->pdg_id[motherpart_idx]);
                   motherpart_idx_tmp = motherpart_idx;
                } 
              
                 mother_pid = mother_pid_tmp;
                 //std::cout << "Final mother ID = " << mother_pid << std::endl;
            }
            
            if( (mother_pid != stats.pmap.at("MuonMotherIDs").first) && (mother_pid != stats.pmap.at("MuonMotherIDs").second) ) continue;
         }
      }

      active_part->at(genMaper.at(particle_id)->ePos)->push_back(j);
      //if(genMaper.at(particle_id)->ePos == CUTS::eGTau) std::cout << "index of Gen collection stored = " << j << std::endl;
    }

  }
}

// --- Function that applies selections to hadronic taus at gen-level (stored in the GenVisTau list) --- //
void Analyzer::getGoodGenHadronicTaus(const PartStats& stats){

  // Loop over all gen-level hadronic taus stored in the corresponding list to apply certain selections
  for(size_t i=0; i < _GenHadTau->size(); i++){

    if(_GenHadTau->pt(i) < stats.pmap.at("HadTauPtCut").first || _GenHadTau->pt(i) > stats.pmap.at("HadTauPtCut").second || abs(_GenHadTau->eta(i)) > stats.dmap.at("HadTauEtaCut")) continue;
    else if( stats.bfind("DiscrByTauDecayMode") && (_GenHadTau->decayMode[i] < stats.pmap.at("TauDecayModes").first || _GenHadTau->decayMode[i] > stats.pmap.at("TauDecayModes").second)) continue;

    active_part->at(CUTS::eGHadTau)->push_back(i);
  }
}

// --- Function that applies selections to jets at gen-level (stored in the GenJets list) --- //
 // The jet flavour is determined based on the "ghost" hadrons clustered inside a jet.
 // More information about Jet Parton Matching at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools

 void Analyzer::getGoodGenJets(const PartStats& stats){

   // Loop over all gen-level jets from the GenJet collection to apply certain selections
   for(size_t i=0; i < _GenJet->size(); i++){
     if(stats.bfind("DiscrJetByPtandEta")){
       if(_GenJet->pt(i) < stats.pmap.at("JetPtCut").first || _GenJet->pt(i) > stats.pmap.at("JetPtCut").second || abs(_GenJet->eta(i)) > stats.dmap.at("JetEtaCut")) continue;
     }
     else if(stats.bfind("DiscrByPartonFlavor")){
       int jetPartonFlavor = _GenJet->genPartonFlavor[i]; // b-jet: jetPartonFlavor = 5, c-jet: jetPartonFlavor = 4, light-jet: jetPartonFlavor = 1,2,3,21, undefined: jetPartonFlavor = 0

       if(abs(jetPartonFlavor) < stats.pmap.at("PartonFlavorRange").first || abs(jetPartonFlavor) > stats.pmap.at("PartonFlavorRange").second) continue;
     }
     active_part->at(CUTS::eGJet)->push_back(i);
   }
 }

 // --- Function that applies selections to b-jets at gen-level (stored in the GenVisTau list) --- //
 void Analyzer::getGoodGenBJets(const PartStats& stats){

   // Loop over all gen-level jets from the GenJet collection to apply certain selections
   for(size_t i=0; i < _GenJet->size(); i++){
    //std::cout << "GenJet #" << i << std::endl;

    int jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[i]); // b-jet: jetFlavor = 5, c-jet: jetFlavor = 4, light-jet: jetFlavor = 0
    //std::cout << "GenJet hadron flavor = " << jetHadronFlavor << ", parton flavor = " << _GenJet->genPartonFlavor[i] << std::endl;
    // Only consider those jets that have parton/hadron flavor = 5
    if(abs(_GenJet->genPartonFlavor[i]) != 5 || abs(jetHadronFlavor) != 5) continue;
    // std::cout << "This is a b-jet!!" << std::endl;

    if(stats.bfind("DiscrBJetByPtandEta")){
      // std::cout << "Applied a gen level cut to the pt and eta of the b-jets" << std::endl;
      // std::cout << "b-jet pt = " << _GenJet->pt(i) << ", eta = " << _GenJet->eta(i) << std::endl;
      if(_GenJet->pt(i) < stats.pmap.at("BJetPtCut").first || _GenJet->pt(i) > stats.pmap.at("BJetPtCut").second || abs(_GenJet->eta(i)) > stats.dmap.at("BJetEtaCut")) continue;
      // std::cout << "Passed kinematic cuts." << std::endl;
    }

    if(stats.bfind("DiscrBJetByMotherID")){

      // std::cout << "Discriminating by mother particle ID" << std::endl;
      TLorentzVector gen_jet = _GenJet->p4(i);
      TLorentzVector genparticle(0,0,0,0);
      int nmatches = 0;

      for(size_t j = 0; j < _Gen->size(); j++){
        // Find the same particle but in the full gen particle collection, to have access to the mother particle ID information.
        genparticle = _Gen->p4(j);
        // std::cout << "Gen particle #" << j << std::endl;
        // std::cout << "Gen particle pt = " << genparticle.Pt() << ", eta = " <<  genparticle.Eta() << ", pdg_id = " << _Gen->pdg_id[j] << ", mother ID = " << _Gen->pdg_id[_Gen->genPartIdxMother[j]] << std::endl;
        // Check through a stringent deltaR requirement
        /*
        if(genparticle.DeltaR(gen_jet) < stats.dmap.at("BJetQuarkDeltaR") && abs(_Gen->pdg_id[j]) == 5 && abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]) == 6){
          std::cout << "This particle matches the original b-jet" << std::endl;
          nmatches++;
        }
        */
        if(genparticle.DeltaR(gen_jet) < stats.dmap.at("BJetQuarkDeltaR") && abs(_Gen->pdg_id[j]) == 5){
          // Get the ID of the mother particle
          int motherPID = _Gen->pdg_id[_Gen->genPartIdxMother[j]];
          // std::cout << "ID of matched gen particle = " << _Gen->pdg_id[j] << ", and its mother particle = " << motherPID << std::endl;

          if( abs(motherPID) == 6){
            // std::cout << "Directly matched to a top quark decay." << std::endl;
            nmatches++;
            break;
          }
          else if ( abs(motherPID) != 6){
            // std::cout << "Not directly matched to a top quark decay. Checking decay chain..." << std::endl;

            int motheridx = _Gen->genPartIdxMother[j];
            int mPID_tmp = motherPID;
            // std::cout << "Initial mother particle ID = " << motherPID << std::endl;
            int decaylevel = 0;
            while(abs(mPID_tmp) != 6){
              decaylevel++;
              // Get the index of the mother particle and assign the new mother PID
              int motheridx_tmp = _Gen->genPartIdxMother[motheridx];
              mPID_tmp = _Gen->pdg_id[motheridx_tmp];
              // std::cout << "mother ID (" << decaylevel << " generations before) = " << mPID_tmp << std::endl;
              motheridx = motheridx_tmp;
              if(abs(mPID_tmp) == 6){
                // std::cout << "Found it's coming from a top quark! " << std::endl;
                nmatches++;
                break;
              }
              else if(mPID_tmp == 0){
                // std::cout << "Never matched to a top quark :( " << std::endl;
                break;
              }
            }
          }
        }
      }
      // std::cout << "nmatches = " << nmatches << std::endl;
      //std::cout << "genparticle mother ID = " << _Gen->pdg_id[_Gen->genPartIdxMother[gen_index]] << ", pdg_id = " << _Gen->pdg_id[gen_index] << std::endl;
      // Get mother particle ID:
      if(nmatches != 1) continue;
      // std::cout << "This b-quark is coming from a top quark! " << std::endl;
    }
    // std::cout << "Gen b-jet vector with index = " << i << " stored in eGBJet." << std::endl;
    active_part->at(CUTS::eGBJet)->push_back(i);
   }
 }

// --- Function that gets the Lorentz vector of taus that decayed hadronically using the tagging method in Analyzer::getGenHadronicTauNeutrinos() --- //
TLorentzVector Analyzer::getGenVisibleTau4Vector(int gentau_idx, int gentaunu_idx){
  TLorentzVector visTau4Vector(0,0,0,0);

  if(gentaunu_idx != -1 && gentau_idx != -1){
    visTau4Vector = _Gen->p4(gentau_idx) - _Gen->p4(gentaunu_idx);
  }
  return visTau4Vector;
}

// --- Function that looks for tau neutrinos coming from hadronic tau decays in the full gen-level particle list --- //
void Analyzer::getGoodGenHadronicTauNeutrinos(const PartStats& stats){
  int genTauNeutrino_idx = -1, genHadTauNeutrino_idx = -1;   // integers for the particle indices in the gen-particle vector.
  TLorentzVector genVisHadTau4Vector(0,0,0,0);

  // Loop over the gen-level taus that satisfied the conditions imposed in getGoodGen, which are stored at CUTS::eGTau.
  for(auto gentau_it : *active_part->at(CUTS::eGTau)){ // (gentau_it) is the index of the gen-level tau in the gen-particles vector.
    // For each iteration, reset the tau neutrino index and the leptonic flag.
    genTauNeutrino_idx = -1, genHadTauNeutrino_idx = -1;
    // isTauLeptonicDecay = false;

    // Make sure that the status code of the current gen-tau is 2, meaning, it's an unstable particle about to decay
    if(_Gen->status[(gentau_it)] != 2) continue;

    // Loop over all gen-level particles to find those that come from the decay of the current gen-tau.
    for(size_t genpart_idx = 0; genpart_idx < _Gen->size(); genpart_idx++) {

      // Check that the mother particle index matches the index of the current gen-tau and that the particle is any neutrino:
      if( abs(_Gen->genPartIdxMother[genpart_idx]) != (gentau_it)) continue; // || abs(_Gen->pdg_id[genpart_idx]) != 12 || abs(_Gen->pdg_id[genpart_idx]) != 14 || abs(_Gen->pdg_id[genpart_idx]) != 16) continue;

      // Look specifically for electron, muon and tau neutrinos. They are enough to tell if a decay was hadronic or leptonic.
      // If in the particle decays there is an electron or muon neutrino, the decay is leptonic and we break this for loop and continue with the next gen-tau in CUTS::eGTau
      if( ( abs(_Gen->pdg_id[genpart_idx]) == 12 || ( abs(_Gen->pdg_id[genpart_idx]) == 14 ) ) ){
        genTauNeutrino_idx = -1;
        break;
      }
      // Since we reject all taus that have e/mu neutrinos, we will only get hadronic tau decays and we want to store the index of the neutrino for these events.
      // As soon as it's found, then break the loop to not go over the entire list of gen particles for this events.
      else if( abs(_Gen->pdg_id[genpart_idx]) == 16){ // genTauNeutrino_idx = genpart_idx;
        genTauNeutrino_idx = genpart_idx;
      }
    }

    // Check if the current tau decayed hadronically or leptonically and assign the index accordingly
    genHadTauNeutrino_idx = genTauNeutrino_idx;

    // Apply kinematic cuts on the visible hadronic tau vector
    genVisHadTau4Vector = getGenVisibleTau4Vector(gentau_it, genHadTauNeutrino_idx);
    if(genVisHadTau4Vector.Pt() < stats.pmap.at("HadTauPtCut").first || genVisHadTau4Vector.Pt() > stats.pmap.at("HadTauPtCut").second || abs(genVisHadTau4Vector.Eta()) > stats.dmap.at("HadTauEtaCut")){
      genHadTauNeutrino_idx = -1;
    }

    active_part->at(CUTS::eGNuTau)->push_back(genHadTauNeutrino_idx);
  }

}
/*
void Analyzer::getGoodGenBJet() { //01.16.19
  for (size_t j=0; j < _Gen->size(); j++){
    int id = abs(_Gen->pdg_id[j]);
    int motherid = abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]);
    int motherind = abs(_Gen->genPartIdxMother[j]);
    if(id == 5 && motherid == 6)
      {for(size_t k=0; k < _Gen->size(); k++){
    if(abs(_Gen->pdg_id[k]) == 24 && _Gen->genPartIdxMother[k] == motherind)
      {active_part->at(CUTS::eGBJet)->push_back(j);
      }
  }
      }
  }
}
*/

double Analyzer::getTopBoostWeight(){ //01.15.19
  double topPt; //initialize a value to hold the p_T of the top
  double topBarPt; //initialize a value to hold the p_T of the tbar
  double SFtop = 1; //initialize a value to hold the top SF
  double SFtopBar = 1; //initialize a value to hold the tbar SF
  double SFttbar = 1; //holds SF for both

  for(size_t j = 0; j < _Gen->size(); j++) {
    int id = _Gen->pdg_id[j]; //grab the particle ID
    int daught = _Gen->numDaught[j]; //grab the number of daughter particles
    if(id == 6 && daught == 2){ //check that a top has two daughters
      topPt = _Gen->pt(j); //grab its p_T
      SFtop = exp(0.0615 - (0.0005 * topPt));} //calculate the top SF
    if(id == -6 && daught == 2){ //check that a tbar has two daughters
      topBarPt = _Gen->pt(j); //grab its p_T
      SFtopBar = exp(0.0615 - (0.0005 * topBarPt));} //calculate the tbar SF
    SFttbar = sqrt(SFtop * SFtopBar);} //calculate the total SF

  return SFttbar;
}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) {
    active_part->at(ePos) = goodParts[ePos];
    return;
  }

  int i = 0;
  for(auto lvec: lep) {

    bool passCuts = true;
    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) passCuts = passCuts && false;
    else if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) passCuts = passCuts && false;

    if((lep.pstats.at("Smear").bfind("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    for( auto cut: stats.bset) {
      if(!passCuts) break;
      else if(cut == "DoDiscrByIsolation") {
        double firstIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
        double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : stats.bfind("FlipIsolationRequirement");
        passCuts = passCuts && lep.get_Iso(i, firstIso, secondIso);
      }
      else if(cut == "DiscrIfIsZdecay" && lep.type != PType::Tau ) passCuts = passCuts && isZdecay(lvec, lep);
      else if(cut == "DiscrByMetDphi") passCuts = passCuts && passCutRange(absnormPhi(lvec.Phi() - _MET->phi()), stats.pmap.at("MetDphiCut"));
      else if(cut == "DiscrByMetMt") passCuts = passCuts && passCutRange(calculateLeptonMetMt(lvec), stats.pmap.at("MetMtCut"));
      /////muon cuts
      else if(lep.type == PType::Muon){
        if(cut == "DoDiscrByTightID") passCuts = passCuts && _Muon->tight[i];
        else if(cut == "DoDiscrBySoftID") passCuts = passCuts && _Muon->soft[i];
        else if(cut == "DoDiscrBySUSYSoftID"){
        	passCuts = passCuts && discrBySUSYSoftMuons(i, stats.bfind("ApplyLooseNotTightID"), stats.bfind("ApplyTightID"));
        }
      }
      ////electron cuts
      else if(lep.type == PType::Electron){
        if(cut == "DoDiscrByCBID"){
          std::bitset<8> idvariable(_Electron->cutBased[i]);
          if(ival(ePos) - ival(CUTS::eRElec2)){ //test if it is electron1
            passCuts = passCuts && (_Electron->cbIDele1&idvariable).count();
          }else{
            passCuts = passCuts && (_Electron->cbIDele2&idvariable).count();
          }
        }
        else if(cut == "DoDiscrBymvaID"){
          if(stats.bfind("DiscrBymvaWP80")) passCuts = passCuts && _Electron->mvaFall17V2Iso_WP80[i];

          else if(stats.bfind("DiscrBymvaWP90")) passCuts = passCuts && _Electron->mvaFall17V2Iso_WP90[i];

          else if(stats.bfind("DiscrBymvaWPL")) passCuts = passCuts && _Electron->mvaFall17V2Iso_WPL[i];
        }
        else if(cut == "DoDiscrByHEEPID"){
        	passCuts = passCuts && _Electron->isPassHEEPId[i];
        }
        else if(cut == "DoDiscrBySUSYSoftID"){
        	passCuts = passCuts && discrBySUSYSoftElectrons(i, stats.bfind("ApplyLooseNotTightID"), stats.bfind("ApplyTightID"));
        }
      }
      else if(lep.type == PType::Tau){
        if(cut == "DoDiscrByCrackCut") passCuts = passCuts && !isInTheCracks(lvec.Eta());
        /////tau cuts
        else if(cut == "DoDzCut") passCuts = passCuts && (abs(_Tau->dz[i]) <= stats.dmap.at("DzCutThreshold"));
        else if(cut == "DoDiscrByLeadTrack") passCuts = passCuts && (_Tau->leadTkPtOverTauPt[i]*_Tau->pt(i) >= stats.dmap.at("LeadTrackThreshold"));
             // ----Electron and Muon vetos
        else if (cut == "DoDiscrAgainstElectron") passCuts = passCuts && _Tau->pass_against_Elec(ePos, i);
        else if (cut == "SelectTausThatAreElectrons") passCuts = passCuts && !_Tau->pass_against_Elec(ePos, i);

        else if (cut == "DoDiscrAgainstMuon") passCuts = passCuts && _Tau->pass_against_Muon(ePos, i);
        else if (cut == "SelectTausThatAreMuons") passCuts = passCuts &&  !_Tau->pass_against_Muon(ePos, i);

        else if(cut == "DiscrByProngType") {
          passCuts = passCuts && (stats.smap.at("ProngType").find("hps") == std::string::npos || _Tau->DecayModeNewDMs[i] != 0);
          // passCuts = passCuts && passProng(stats.smap.at("ProngType"), _Tau->decayMode[i]); //original.
          passCuts = passCuts && passProng(stats.smap.at("ProngType"), _Tau->decayModeInt[i]);
        }
        else if(cut == "decayModeFindingNewDMs") passCuts = passCuts && _Tau->DecayModeNewDMs[i] != 0;
        // else if(cut == "decayModeFinding") passCuts = passCuts && _Tau->DecayMode[i] != 0; // original
        else if(cut == "decayModeFinding") passCuts = passCuts && _Tau->DecayModeOldDMs[i] != 0;
              // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      }
      else std::cout << "cut: " << cut << " not listed" << std::endl;
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }

  return;
}

bool Analyzer::discrBySUSYSoftElectrons(int index, bool applyloosenottight, bool applytight){

	bool passCuts = true;

	// std::cout << "Electron pt = " << _Electron->pt(index) << ", |eta| = " << fabs(_Electron->eta(index)) << ", dxy = " << _Electron->dxy[index] << ", dz = " << _Electron->dz[index] << std::endl;
	// Apply kinematic cuts (pt and eta)
	passCuts = passCuts && _Electron->pt(index) >= 5.0;
	passCuts = passCuts && fabs( _Electron-> eta(index) ) <= 2.5;

	passCuts = passCuts && fabs(_Electron->dxy[index]) < 0.05;
	passCuts = passCuts && fabs(_Electron->dz[index]) < 0.1;

	// std::cout << "passCuts = " << passCuts << std::endl;

	if(applyloosenottight){

		// std::cout << "applying loosenottight" << std::endl;
		// std::cout << "Electron ip3d = " << _Electron->ip3d[index] << ", sip3d = " << _Electron->sip3d[index] << ", abs. iso. = " << (_Electron->pfRelIso03_all[index] * _Electron->pt(index)) << " (req < " << (20.0 + (300.0 / _Electron->pt(index) ) ) << "), rel. iso = " << _Electron->pfRelIso03_all[index] << std::endl;

		// Selection: 3D impact parameter
		passCuts = passCuts && _Electron->ip3d[index] < 0.0175;
		passCuts = passCuts && _Electron->sip3d[index] < 2.5;

		// Isolation requirements
		passCuts = passCuts && ((_Electron->pfRelIso03_all[index] * _Electron->pt(index)) < (20.0 + (300.0 / _Electron->pt(index) ) ) ); // absolute isolation
		passCuts = passCuts && _Electron->pfRelIso03_all[index] < 1.0; // relative isolation

		// Electron MVA custom ID
		// std::cout << "Applying custom electron MVA ID "<< std::endl << _Electron->customSoftLooseEleMVAId(index, _Electron->p4(index), runyear) << std::endl;
		passCuts = passCuts && _Electron->customSoftLooseEleMVAId(index, _Electron->p4(index), runyear);

		// std::cout << "passCuts = " << passCuts << std::endl;
	}
	else if(applytight){

		// std::cout << "applying tight" << std::endl;
		// std::cout << "Electron ip3d = " << _Electron->ip3d[index] << ", sip3d = " << _Electron->sip3d[index] << ", abs. iso. = " << (_Electron->pfRelIso03_all[index] * _Electron->pt(index)) << " (req < " << (20.0 + (300.0 / _Electron->pt(index) ) ) << "), rel. iso = " << _Electron->pfRelIso03_all[index] << std::endl;

		// Selection: 3D impact parameter
		passCuts = passCuts && _Electron->ip3d[index] < 0.01;
		passCuts = passCuts && _Electron->sip3d[index] < 2.0;

		// Isolation requirements
		passCuts = passCuts && ((_Electron->pfRelIso03_all[index] * _Electron->pt(index)) < 5.0 ); // absolute isolation
		passCuts = passCuts && _Electron->pfRelIso03_all[index] < 0.5; // relative isolation

		// Electron MVA custom ID
		// std::cout << "Applying custom electron MVA ID "<< std::endl << _Electron->customSoftTightEleMVAId(index, _Electron->p4(index), runyear) << std::endl;
		passCuts = passCuts && _Electron->customSoftTightEleMVAId(index, _Electron->p4(index), runyear);

		// Additional b-jet veto.
		int elejetidx = _Electron->associatedJetIndex[index]; // get the index of the jet collection which is associated with the current jet.
		// std::cout << "Associated jet index = " << elejetidx << std::endl;
		// Apply a veto to this jet to which we further apply pt and eta cuts and b-tagging cuts (Loose WPs from DeepCSV in all years, custom in 2016):

		if(elejetidx != -1){
			// std::cout << "Associated jet pt = " << _Jet->pt(elejetidx) << ", |eta| = " << fabs(_Jet->eta(elejetidx)) << ", b-tag discr. = " << _Jet->bDiscriminatorDeepCSV[elejetidx] << std::endl;
			bool associatedBJet = _Jet->pt(elejetidx) > _Jet->pstats["BJet"].dmap.at("PtCut") && passCutRange(fabs(_Jet->eta(elejetidx)), _Jet->pstats["BJet"].pmap.at("EtaCut")) && _Jet->bDiscriminatorDeepCSV[elejetidx] > susysoftlepidloosebtagdiscr.at(runyear);
			// std::cout << "Is this an associated b-jet? " << associatedBJet << std::endl;
			passCuts = passCuts && !associatedBJet;
		}

		// std::cout << "passCuts = " << passCuts << std::endl;
	}

	// Missing pixel hits
	int missingPixelHits = static_cast<unsigned>(_Electron->lostHits[index]);
	// std::cout << "Electron lost hits = " << missingPixelHits << ", conversion veto = " << _Electron->conversionVeto[index] << std::endl;

	passCuts = passCuts && (missingPixelHits == 0);
	// Conversion vertex veto
	passCuts = passCuts && _Electron->conversionVeto[index];

	// std::cout << "passCuts = " << passCuts << std::endl;
	// std::cout << "------------" << std::endl;
	return passCuts;
}

bool Analyzer::discrBySUSYSoftMuons(int index, bool applyloosenottight, bool applytight){

	bool passCuts = true;
	// Apply kinematic cuts (pt and eta)
	// std::cout << "Muon pt = " << _Muon->pt(index) << ", |eta| = " << fabs(_Muon->eta(index)) << ", dxy = " << _Muon->dxy[index] << ", dz = " << _Muon->dz[index] << std::endl;
	passCuts = passCuts && _Muon->pt(index) >= 3.5;
	passCuts = passCuts && fabs( _Muon-> eta(index) ) <= 2.4;

	passCuts = passCuts && fabs(_Muon->dxy[index]) < 0.05;
	passCuts = passCuts && fabs(_Muon->dz[index]) < 0.1;

	// std::cout << "passCuts = " << passCuts << std::endl;

	if(applyloosenottight){
		// std::cout << "applying loosenottight" << std::endl;

		// std::cout << "Muon ip3d = " << _Muon->ip3d[index] << ", sip3d = " << _Muon->sip3d[index] << ", abs. iso. = " << (_Muon->pfRelIso03_all[index] * _Muon->pt(index)) << " (req < " << (20.0 + (300.0 / _Muon->pt(index) ) ) << "), rel. iso = " << _Muon->pfRelIso03_all[index] << std::endl;
		// Selection: 3D impact parameter
		passCuts = passCuts && _Muon->ip3d[index] < 0.0175;
		passCuts = passCuts && _Muon->sip3d[index] < 2.5;

		// Isolation requirements
		passCuts = passCuts && ((_Muon->pfRelIso03_all[index] * _Muon->pt(index)) < (20.0 + (300.0 / _Muon->pt(index) ) ) ); // absolute isolation
		passCuts = passCuts && _Muon->pfRelIso03_all[index] < 1.0; // relative isolation

		// std::cout << "passCuts = " << passCuts << std::endl;
	}
	else if(applytight){

		// std::cout << "applying tight" << std::endl;
		// std::cout << "Muon ip3d = " << _Muon->ip3d[index] << ", sip3d = " << _Muon->sip3d[index] << ", abs. iso. = " << (_Muon->pfRelIso03_all[index] * _Muon->pt(index)) << "(req < 5), rel. iso = " << _Muon->pfRelIso03_all[index] << std::endl;

		// Selection: 3D impact parameter
		passCuts = passCuts && _Muon->ip3d[index] < 0.01;
		passCuts = passCuts && _Muon->sip3d[index] < 2.0;

		// Isolation requirements
		passCuts = passCuts && ((_Muon->pfRelIso03_all[index] * _Muon->pt(index)) < 5.0 ); // absolute isolation
		passCuts = passCuts && _Muon->pfRelIso03_all[index] < 0.5; // relative isolation

		// Additional b-jet veto.
		int mujetidx = _Muon->associatedJetIndex[index]; // get the index of the jet collection which is associated with the current jet.
		// std::cout << "Associated jet index = " << mujetidx << std::endl;
		if(mujetidx != -1){
			// std::cout << "Associated jet pt = " << _Jet->pt(mujetidx) << ", |eta| = " << fabs(_Jet->eta(mujetidx)) << ", b-tag discr. = " << _Jet->bDiscriminatorDeepCSV[mujetidx] << std::endl;
			bool associatedBJet = _Jet->pt(mujetidx) > _Jet->pstats["BJet"].dmap.at("PtCut") && passCutRange(fabs(_Jet->eta(mujetidx)), _Jet->pstats["BJet"].pmap.at("EtaCut")) && _Jet->bDiscriminatorDeepCSV[mujetidx] > susysoftlepidloosebtagdiscr.at(runyear);
			// std::cout << "Is this an associated b-jet? " << associatedBJet << std::endl;
			passCuts = passCuts && !associatedBJet;
		}
		// std::cout << "passCuts = " << passCuts << std::endl;
	}

	// Soft and loose ID muon POG requirements
	// std::cout << "Additional ID requirements: passLoose = " << _Muon->looseId[index] << ", passSoft = " << _Muon->soft[index] << std::endl;
	passCuts = passCuts && _Muon->looseId[index];
	passCuts = passCuts && _Muon->soft[index];

	// std::cout << "passCuts = " << passCuts << std::endl;
	// std::cout << "----------" << std::endl;

	return passCuts;
}



bool Analyzer::passJetVetoEEnoise2017(int jet_index){
   // Get the jet pT and raw factors to calculate the raw jet pt:
   TLorentzVector jet_RecoP4 = _Jet->RecoP4(jet_index);
   double jet_rawFactor = _Jet->rawFactor[jet_index];

   double jet_rawPt = jet_RecoP4.Pt() * (1.0 - jet_rawFactor);

   //std::cout << "Jet #" << jet_index << ", raw pt = " << jet_rawPt << ", eta = " << _Jet->eta(jet_index) << ", reco eta = " << jet_RecoP4.Eta() << std::endl;

   // Check if this jet is in the problematic pt-eta region: if yes, then jet does not pass the veto requirement
   if(jet_rawPt < 50.0 && (abs(jet_RecoP4.Eta()) > 2.65 && abs(jet_RecoP4.Eta()) < 3.139)){
    //std::cout << "Jet in the EE noisy region, throwing it out." << std::endl;
    return false;
   }
   // Otherwise, return true (passes the veto)
   //std::cout << "Jet not in the EE noise region, keeping it." << std::endl;
   return true;
 }

 void Analyzer::getAdditionalBadEENoiseRecoJets(const PartStats& stats) {

   int i=0;

   for(auto lvec: *_Jet) {

     bool passCuts = true;
     double dphi1rjets = normPhi(lvec.Phi() - _MET->phi());

     passCuts = passCuts && ( (fabs(lvec.Eta()) > 2.65) && (fabs(lvec.Eta()) < 3.14) ); // Eta region for EE noise
     passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut") && lvec.Pt() < 80.0);

     for( auto cut: stats.bset) {
       if(!passCuts) break;
       else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(i);
       else if(cut == "ApplyLooseID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut"));
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedLooseJetID(i);
         }
       }
       else if(cut == "ApplyTightID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedTightJetID(i);
         }
       }
       else if(cut == "ApplyTightLepVetoID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedTightLepVetoJetID(i);
         }
       }
       else if(cut == "RemoveOverlapWithJs") passCuts = passCuts && !isOverlapingC(lvec, *_FatJet, CUTS::eRWjet, stats.dmap.at("JMatchingDeltaR"));
       else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(lvec, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
       // ----anti-overlap requirements
       else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
       else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
       else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
       else if(cut == "DiscrByDPhiMet") passCuts = passCuts && passCutRange(fabs(dphi1rjets), stats.pmap.at("DPhiMetCut")); //01.17.19
     }

     if(passCuts) additionalEEnoiseJets.push_back(i);
     i++;

   }
 }

 void Analyzer::passAdditionalEENoiseVeto(CUTS ePos, const PartStats& stats, std::string year) {

  if(! neededCuts.isPresent(ePos)) return;

  if(year.compare("2017") != 0 || distats["Run"].bfind("ApplyAdditionalEEnoiseVeto2017") == 0){
    active_part->at(ePos)->push_back(0);
    return;
  }

   int i=0;

   std::vector<int> jetsminuseta3p2to2p6;
   std::vector<int> jetspluseta2p6to3p2;
   std::vector<int> jetsminuseta4p7to0;
   std::vector<int> jetspluseta0to4p7;

   for(auto lvec: *_Jet) {

     bool passCuts = true;

     passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut"));

     for( auto cut: stats.bset) {
       if(!passCuts) break;
       else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(i);
       else if(cut == "ApplyLooseID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut"));
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedLooseJetID(i);
         }
       }
       else if(cut == "ApplyTightID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedTightJetID(i);
         }
       }
       else if(cut == "ApplyTightLepVetoID"){
         if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
           if(!stats.bfind("FailPUJetID")){
             passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
           } else {
             passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
           }
         } else {
           passCuts = passCuts && _Jet->passedTightLepVetoJetID(i);
         }
       }
       else if(cut == "RemoveOverlapWithJs") passCuts = passCuts && !isOverlapingC(lvec, *_FatJet, CUTS::eRWjet, stats.dmap.at("JMatchingDeltaR"));
       else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(lvec, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
       // ----anti-overlap requirements
       else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
       else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
       else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
       else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
     }

     if(passCuts){
      // Check if it is in the positive or negative eta region
      if(lvec.Eta() < 0){ jetsminuseta4p7to0.push_back(i);  }
      else if(lvec.Eta() > 0){ jetspluseta0to4p7.push_back(i); }

      if( (lvec.Eta() > -3.15) && (lvec.Eta() < -2.66)){ jetsminuseta3p2to2p6.push_back(i); }
      if( (lvec.Eta() > 2.66) && (lvec.Eta() > 3.15 )){ jetspluseta2p6to3p2.push_back(i); }
     }

     i++;

   }

   bool passVeto = true;
   // Check the different cases possible for this veto
   if( (jetsminuseta3p2to2p6.size() > 0) && (jetspluseta2p6to3p2.size() == 0)){
      if(jetspluseta0to4p7.size() == 0){
        passVeto = false;
      }
   } else if( (jetsminuseta3p2to2p6.size() == 0) && (jetspluseta2p6to3p2.size() > 0)){
      if(jetsminuseta4p7to0.size() == 0){
        passVeto = false;
      }
   } else if( (jetsminuseta3p2to2p6.size() == 1) && (jetspluseta2p6to3p2.size() == 1)){
    // Apply a deltaPt cut
    TLorentzVector deltaP = _Jet->p4(jetspluseta2p6to3p2.at(0)) - _Jet->p4(jetsminuseta3p2to2p6.at(0));

    float ratio_deltaPtHT = deltaP.Pt() / _MET->HT();

    if(ratio_deltaPtHT < 0.4) passVeto = false;
   }

   if(passVeto){
      active_part->at(ePos)->push_back(0);
      return;
   }

}


////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats, const int syst) {

  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_Jet) {

    bool passCuts = true;
    double dphi1rjets = normPhi(lvec.Phi() - _MET->phi());
    if(ePos == CUTS::eRCenJet) passCuts = passCuts && (fabs(lvec.Eta()) < 2.5);
    else  passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    for( auto cut: stats.bset) {
      if(!passCuts) break;
      else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(i);
    	/// BJet specific
      else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyLooseID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
            passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut"));
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        } else {
          passCuts = passCuts && _Jet->passedLooseJetID(i);
        }
      }
      else if(cut == "ApplyTightID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
            passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        } else {
          passCuts = passCuts && _Jet->passedTightJetID(i);
        }
      }
      else if(cut == "ApplyTightLepVetoID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
            passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        } else {
          passCuts = passCuts && _Jet->passedTightLepVetoJetID(i);
        }
      }
      else if(cut == "RemoveOverlapWithJs") passCuts = passCuts && !isOverlapingC(lvec, *_FatJet, CUTS::eRWjet, stats.dmap.at("JMatchingDeltaR"));
      else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(lvec, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
      // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
      else if(cut == "DiscrByDPhiMet") passCuts = passCuts && passCutRange(fabs(dphi1rjets), stats.pmap.at("DPhiMetCut")); //01.17.19
    }
    if(_Jet->pstats["BJet"].bfind("RemoveBJetsFromJets") and ePos!=CUTS::eRBJet){
      passCuts = passCuts && find(active_part->at(CUTS::eRBJet)->begin(), active_part->at(CUTS::eRBJet)->end(), i) == active_part->at(CUTS::eRBJet)->end();
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;

  }
}

void Analyzer::getGoodRecoLeadJets(CUTS ePos, const PartStats& stats, const int syst) {

  if(!neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  // Clean jetPtIndexVector after looping over first and second leading jets for each iteration.
  if(ePos == CUTS::eR1stJet && jetPtIndexVector.size() != 0){
      jetPtIndexVector.clear();
      jetPtIndexVector.shrink_to_fit();
  }

  //note the leading jet has to be selected first!
  //Do this only once for the first leading jet:
  if(ePos == CUTS::eR1stJet){
    for(auto it : *active_part->at(CUTS::eRJet1)) {
      jetPtIndexVector.push_back(std::make_pair(_Jet->pt(it),it));
    }
    sort(jetPtIndexVector.begin(),jetPtIndexVector.end());
  }

  if(ePos == CUTS::eR1stJet && jetPtIndexVector.size()>0){

    bool passCuts = true;
    int i = jetPtIndexVector.at(jetPtIndexVector.size()-1).second;

    if(stats.bfind("DoDiscrByThisJet")){

      TLorentzVector leadjetp4 = _Jet->p4(i);
      double dphi1rjets = normPhi(leadjetp4.Phi() - _MET->phi());
      // Applying cuts for FirstLeadingJet block
      passCuts = passCuts && passCutRange(fabs(leadjetp4.Eta()), stats.pmap.at("EtaCut"));
      passCuts = passCuts && (leadjetp4.Pt() > stats.dmap.at("PtCut"));

      for( auto cut: stats.bset) {
        if(!passCuts) break;
        else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(i);
        /// BJet specific
        else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[i] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[i] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[i] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyLooseID"){
          if(stats.bfind("ApplyPileupJetID") && leadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut"));
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedLooseJetID(i);
          }
        }
        else if(cut == "ApplyTightID"){
          if(stats.bfind("ApplyPileupJetID") && leadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedTightJetID(i);
          }
        }
        else if(cut == "ApplyTightLepVetoID"){
          if(stats.bfind("ApplyPileupJetID") && leadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedTightLepVetoJetID(i);
          }
        }
        else if(cut == "RemoveOverlapWithJs") passCuts = passCuts && !isOverlapingC(leadjetp4, *_FatJet, CUTS::eRWjet, stats.dmap.at("JMatchingDeltaR"));
        else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(leadjetp4, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
        // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
        else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
        else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(leadjetp4, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
        else if(cut == "DiscrByDPhiMet") passCuts = passCuts && passCutRange(fabs(dphi1rjets), stats.pmap.at("DPhiMetCut")); //01.17.19
      }
    }

    if(passCuts){
      active_part->at(ePos)->push_back(jetPtIndexVector.back().second);
    }
  }
  else if(ePos == CUTS::eR2ndJet && jetPtIndexVector.size()>1){

    bool passCuts = true;
    int idx = jetPtIndexVector.size() - 2;
    int j = jetPtIndexVector.at(idx).second;

    if(stats.bfind("DoDiscrByThisJet")){

      TLorentzVector subleadjetp4 = _Jet->p4(j);
      double dphi1rjets = normPhi(subleadjetp4.Phi() - _MET->phi());

      // Applying cuts for SecondLeadingJet block
      passCuts = passCuts && passCutRange(fabs(subleadjetp4.Eta()), stats.pmap.at("EtaCut"));
      passCuts = passCuts && (subleadjetp4.Pt() > stats.dmap.at("PtCut"));

      for( auto cut: stats.bset) {
        if(!passCuts) break;
        else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(j);
        /// BJet specific
        else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[j] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[j] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[j] > stats.dmap.at("JetBTaggingCut"));
        else if(cut == "ApplyLooseID"){
          if(stats.bfind("ApplyPileupJetID") && subleadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(j, stats.dmap.at("PUJetIDCut"));
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(j,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedLooseJetID(j);
          }
        }
        else if(cut == "ApplyTightID"){
          if(stats.bfind("ApplyPileupJetID") && subleadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(j, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(j,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedTightJetID(j);
          }
        }
        else if(cut == "ApplyTightLepVetoID"){
          if(stats.bfind("ApplyPileupJetID") && subleadjetp4.Pt() <= 50.0){
            if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(j, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
            } else {
              passCuts = passCuts && (_Jet->getPileupJetID(j,0) == 0);
            }
          } else {
            passCuts = passCuts && _Jet->passedTightLepVetoJetID(j);
          }
        }
        else if(cut == "RemoveOverlapWithJs") passCuts = passCuts && !isOverlapingC(subleadjetp4, *_FatJet, CUTS::eRWjet, stats.dmap.at("JMatchingDeltaR"));
        else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(subleadjetp4, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
        // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
        else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
        else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(subleadjetp4, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
        else if(cut == "DiscrByDPhiMet") passCuts = passCuts && passCutRange(fabs(dphi1rjets), stats.pmap.at("DPhiMetCut")); //01.17.19
      }
    }

    if(passCuts){
      active_part->at(ePos)->push_back(jetPtIndexVector.at(idx).second);
    }
  }
}


void Analyzer::getGoodRecoBJets(CUTS ePos, const CUTS eGenPos, const PartStats& stats, const int syst) { //01.16.19
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;
  for(auto lvec: *_Jet) {
    bool passCuts = true;
    if( ePos == CUTS::eRCenJet) passCuts = passCuts && (fabs(lvec.Eta()) < 2.5);
    else  passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    for( auto cut: stats.bset) {
      if(!passCuts) break;
      else if(cut == "Apply2017EEnoiseVeto") passCuts = passCuts && passJetVetoEEnoise2017(i);
      /// BJet specific
      else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "MatchBToGen" && !isData){
        int matchedGenJetIndex = _Jet->genJetIdx[i];
        passCuts = passCuts && ( std::find(active_part->at(eGenPos)->begin(), active_part->at(eGenPos)->end(), matchedGenJetIndex) != active_part->at(eGenPos)->end() );
      }
      else if(cut == "ApplyLooseID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
            passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut"));
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        }
        else{
          passCuts = passCuts && _Jet->passedLooseJetID(i);
        }
      }
      else if(cut == "ApplyTightID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        } else {
          passCuts = passCuts && _Jet->passedTightJetID(i);
        }
      }
      else if(cut == "ApplyTightLepVetoID"){
        if(stats.bfind("ApplyPileupJetID") && lvec.Pt() <= 50.0){
          if(!stats.bfind("FailPUJetID")){
              passCuts = passCuts && _Jet->getPileupJetID(i, stats.dmap.at("PUJetIDCut")); // Only apply this cut to low pt jets
          } else {
            passCuts = passCuts && (_Jet->getPileupJetID(i,0) == 0);
          }
        } else {
          passCuts = passCuts && _Jet->passedTightLepVetoJetID(i);
        }
      }
      // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

// The function below sets up the information from the right CSV file in the Pileup folder
// to obtain the functions needed to apply b-tagging SF in an automatic way.
void Analyzer::setupBJetSFInfo(const PartStats& stats, std::string year){

   // Always check that the filenames are up to date!!
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
   static std::map<std::string, std::string> csvv2btagsfcsvfiles = {
    {"2016", "CSVv2_Moriond17_B_H.csv"},
    {"2017", "CSVv2_94XSF_WP_V2_B_F.csv"},
    {"2018", "none.csv"} // No CSVv2 available for 2018
   };

   static std::map<std::string, std::string> deepcsvbtagsfcsvfiles = {
    {"2016", "DeepCSV_2016LegacySF_WP_V1.csv"},
    {"2017", "DeepCSV_94XSF_WP_V4_B_F.csv"},
    {"2018", "DeepCSV_102XSF_WP_V1.csv"}
   };

  static std::map<std::string, std::string> deepflavbtagsfcsvfiles = {
    {"2016", "DeepJet_2016LegacySF_WP_V1.csv"},
    {"2017", "DeepFlavour_94XSF_WP_V3_B_F.csv"},
    {"2018", "DeepJet_102XSF_WP_V1.csv"}
  };

   // This will be the default.
  std::string btagalgoname = "DeepCSV";
  std::string btagsffilename = deepcsvbtagsfcsvfiles.begin()->second;

  // std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  // std::cout << "Setting up the b-jet scale factors... " << std::endl;
  // Get the b-tagging algorithm to use to initialize the appropriate csv file:
  if(stats.bfind("ApplyJetBTaggingCSVv2")){
    btagalgoname = "CSVv2";
    btagsffilename = csvv2btagsfcsvfiles[year];
  }
  else if(stats.bfind("ApplyJetBTaggingDeepCSV")){
    btagalgoname = "DeepCSV";
    btagsffilename = deepcsvbtagsfcsvfiles[year];

  }
  else if(stats.bfind("ApplyJetBTaggingDeepFlav")){
    btagalgoname = "DeepFlav";
    btagsffilename = deepflavbtagsfcsvfiles[year];
  }

  try{
    btagcalib = BTagCalibration(btagalgoname, (PUSPACE+"BJetDatabase/"+btagsffilename).c_str());

    if(stats.bfind("UseBtagSF")){
      std::cout << " ---------------------------------------------------------------------- " << std::endl;
      std::cout << "Setting up the b-jet scale factors... " << std::endl;
      std::cout << "B-jet ID algorithm selected: " << btagalgoname << std::endl;
      std::cout << "Selected working point: " << stats.smap.at("JetBTaggingWP") << std::endl;
      std::cout << "B-tagging SF csv file: " << btagsffilename << std::endl;
      std::cout << " ---------------------------------------------------------------------- " << std::endl;
    }
  }
  catch(std::exception& e){
    std::cerr << "ERROR in setupBJetSFInfo: " << e.what() << std::endl;
    std::cout << "\t Setting dummy b-tagging from " << (PUSPACE+"btagging.csv").c_str() << std::endl;
    btagcalib = BTagCalibration("DeepCSV", (PUSPACE+"btagging.csv").c_str());
  }

 // Check BTagCalibrationStandalone.h for more information

   static std::map<std::string, int> bjetflavors = {
    {"bjet", 0}, {"cjet", 1}, {"lightjet", 2},
   };

   bjetflavor = (BTagEntry::JetFlavor) bjetflavors["bjet"];

   static std::map<std::string, int> btagoperatingpoints = {
    {"loose", 0}, {"medium", 1}, {"tight", 2}, {"reshaping", 3}
   };

   b_workingpoint = (BTagEntry::OperatingPoint) btagoperatingpoints[stats.smap.at("JetBTaggingWP").c_str()];

}

double Analyzer::getBJetSF(CUTS ePos, const PartStats& stats) {
  double bjetSFall = 1.0, bjetSFtemp = 1.0;

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  // Load the info from the btaggin reader
  btagsfreader = BTagCalibrationReader(b_workingpoint, "central");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  int matchedGenJetIndex = -1;
  int jetPartonFlavor = 0;
  int jetHadronFlavor = 0;

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{

    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if( active_part->at(CUTS::eRBJet)->size() == 1) {

      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs(jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){

      // Get the SF for the first b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs(jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }

      // Now multiply by the SF of the second b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(1)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs( jetHadronFlavor) == 5){
        bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
      }
      else{
        bjetSFtemp = bjetSFtemp * 1.0;
      }
    }
  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

double Analyzer::getBJetSFResUp(CUTS ePos, const PartStats& stats) {
  double bjetSFall = 1.0, bjetSFtemp = 1.0;

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  // Load the info from the btaggin reader
  btagsfreader = BTagCalibrationReader(b_workingpoint, "up");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  int matchedGenJetIndex = -1;
  int jetPartonFlavor = 0;
  int jetHadronFlavor = 0;

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{

    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 1){

      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]) ;

      if(jetPartonFlavor == 5 || abs( jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){
      // Get the SF for the first b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs( jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }

      // Now multiply by the SF of the second b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(1)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs( jetHadronFlavor) == 5){
        bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
      }
      else{
        bjetSFtemp = bjetSFtemp * 1.0;
      }
    }
  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

double Analyzer::getBJetSFResDown(CUTS ePos, const PartStats& stats) {
  double bjetSFall = 1.0, bjetSFtemp = 1.0;

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  btagsfreader = BTagCalibrationReader(b_workingpoint, "down");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  int matchedGenJetIndex = -1;
  int jetPartonFlavor = 0;
  int jetHadronFlavor = 0;

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{
    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 1){

      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs(jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }

    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){
      // Get the SF for the first b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(0)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs(jetHadronFlavor) == 5){
        bjetSFtemp = btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      }
      else{
        bjetSFtemp = 1.0;
      }

      // Now multiply by the SF of the second b-jet
      matchedGenJetIndex = _Jet->genJetIdx[active_part->at(CUTS::eRBJet)->at(1)];
      jetPartonFlavor = abs(_GenJet->genPartonFlavor[matchedGenJetIndex]);
      jetHadronFlavor = static_cast<unsigned>(_GenJet->genHadronFlavor[matchedGenJetIndex]);

      if(jetPartonFlavor == 5 || abs(jetHadronFlavor) == 5){
        bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
      }
      else{
        bjetSFtemp = bjetSFtemp * 1.0;
      }

    }
  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

////FatJet specific function for finding the number of V-jets that pass the cuts.
void Analyzer::getGoodRecoFatJets(CUTS ePos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_FatJet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_FatJet) {
    bool passCuts = true;
    passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    ///if else loop for central jet requirements
    for( auto cut: stats.bset) {
      if(!passCuts) break;

      else if(cut == "ApplyJetWTagging") passCuts = passCuts && (passCutRange(_FatJet->tau2[i]/_FatJet->tau1[i], stats.pmap.at("JetTau2Tau1Ratio")) &&
                                                     passCutRange(_FatJet->PrunedMass[i], stats.pmap.at("JetWmassCut")));
    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));

    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
  }
  return false;
}

bool Analyzer::isOverlapingB(const TLorentzVector& lvec, Jet& overlapper, CUTS ePos, double MatchingDeltaR) { //01.17.19
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
  }
  return false;
}


bool Analyzer::isOverlapingC(const TLorentzVector& lvec, FatJet& overlapper, CUTS ePos, double MatchingDeltaR) {
	for(auto it : *active_part->at(ePos)) {
		if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
	}
	return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(std::string prong, int value) {
  return ( (prong.find("1") != std::string::npos &&  (value<5)) ||
  (prong.find("2") != std::string::npos &&  (value>=5 && value<10)) ||
  (prong.find("3") != std::string::npos && (value>=10 && value<12)) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
  (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
  (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
  (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
  (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(CUTS ePos) {

	if(! neededCuts.isPresent(ePos)) return;

  if(ePos == CUTS::eRTrig1){
  	// Loop over all elements of the trigger decisions vector
  	for(bool decision : trigger1namedecisions){
  		if(decision){
  			// If one element is true (1), then store back the event in the triggers vector
  			active_part->at(ePos)->push_back(0);
  			// Clean up the trigger decisions vector to reduce memory usage and have an empty vector for the next event
  			trigger1namedecisions.clear();
        trigger1namedecisions.shrink_to_fit();
  			// End of the function
  			return;
  		}
  	}
  	// If all the elements of the trigger decisions vector are false, then just clean up the trigger decisions vector to reduce memory usage.
  	trigger1namedecisions.clear();
    trigger1namedecisions.shrink_to_fit();
  }

  if(ePos == CUTS::eRTrig2){
    // Loop over all elements of the trigger decisions vector
    for(bool decision : trigger2namedecisions){
      if(decision){
        // If one element is true (1), then store back the event in the triggers vector
        active_part->at(ePos)->push_back(0);
        // Clean up the trigger decisions vector to reduce memory usage and have an empty vector for the next event
        trigger2namedecisions.clear();
        trigger2namedecisions.shrink_to_fit();
        // End of the function
        return;
      }
    }
    // If all the elements of the trigger decisions vector are false, then just clean up the trigger decisions vector to reduce memory usage.
    trigger2namedecisions.clear();
    trigger2namedecisions.shrink_to_fit();
  }


}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eSusyCom)) return;
  std::string systname = syst_names.at(syst);


  if(systname!="orig"){
    //only jet stuff is affected
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eSusyCom)=goodParts[CUTS::eSusyCom];
      return;
    }
  }

  if(active_part->at(CUTS::eR1stJet)->size()==0 || active_part->at(CUTS::eR2ndJet)->size()==0) return;

  TLorentzVector ljet1 = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
  TLorentzVector ljet2 = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));
  TLorentzVector dijet = ljet1 + ljet2;
  double dphi1 = normPhi(ljet1.Phi() - _MET->phi());
  double dphi2 = normPhi(ljet2.Phi() - _MET->phi());

  bool passCuts = true;
  for(auto cut: stats.bset) {
    if(!passCuts) break;
    else if(cut == "DiscrByMass") passCuts = passCuts && passCutRange(dijet.M(), stats.pmap.at("MassCut"));
    else if(cut == "DiscrByPt") passCuts = passCuts && passCutRange(dijet.Pt(), stats.pmap.at("PtCut"));
    else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(ljet1.Eta() - ljet2.Eta()), stats.pmap.at("DeltaEtaCut"));
    else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(ljet1.Phi() - ljet2.Phi()), stats.pmap.at("DeltaPhiCut"));
    else if(cut == "DiscrByOSEta") passCuts = passCuts && (ljet1.Eta() * ljet2.Eta() < 0);
    else if(cut == "DiscrByR1") passCuts = passCuts && passCutRange(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0)), stats.pmap.at("R1Cut"));
    else if(cut == "DiscrByR2") passCuts = passCuts && passCutRange(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), stats.pmap.at("R2Cut"));
    else if(cut == "DiscrByAlpha") {
      double alpha = (dijet.M() > 0) ? ljet2.Pt() / dijet.M() : -1;
      passCuts = passCuts && passCutRange(alpha, stats.pmap.at("AlphaCut"));
    }
    else if(cut == "DiscrByDphi1") passCuts = passCuts && passCutRange(abs(dphi1), stats.pmap.at("Dphi1Cut"));
    else if(cut == "DiscrByDphi2") passCuts = passCuts && passCutRange(abs(dphi2), stats.pmap.at("Dphi2Cut"));

    else std::cout << "cut: " << cut << " not listed" << std::endl;
  }

  if(passCuts)  active_part->at(CUTS::eSusyCom)->push_back(0);
  return;
}

bool Analyzer::passCutRange(double value, const std::pair<double, double>& cuts) {
  return (value > cuts.first && value < cuts.second);
}

bool Analyzer::passCutRangeAbs(double value, const std::pair<double, double>& cuts) {
  return ((value > cuts.first && value < cuts.second) || (value > (cuts.second * -1) && value < (cuts.first * -1)));
}

//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + _MET->px();
  double py = Tobj.Py() + _MET->py();
  double et = Tobj.Et() + _MET->energy();
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}


/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;

  if(howCalc == "InvariantMass") {
    return (Tobj1 + Tobj2).M();
  }


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*_MET->px() - Tobj2.Px()*_MET->py())/denominator;
    double x2 = (Tobj1.Px()*_MET->py() - Tobj1.Py()*_MET->px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  }

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + _MET->p4()).M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + _MET->px())) - (Tobj2.Px() * (Tobj1.Py() + _MET->py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + _MET->py())) - (Tobj1.Py() * (Tobj2.Px() + _MET->px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}


/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);

  if(!lep1.needSyst(syst) && !lep2.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2, llep;

  for(auto i1 : *active_part->at(ePos1)) {
    part1 = lep1.p4(i1);
    // New: Get the rest mass in GeV according to the lepton type, Jun 26, 2020
    float mass1 = leptonmasses.at(lep1.type);

    for(auto i2 : *active_part->at(ePos2)) {
      if(sameParticle && i2 <= i1) continue;
      part2 = lep2.p4(i2);
      // New: Get the rest mass in GeV according to the lepton type, Jun 26, 2020
      float mass2 = leptonmasses.at(lep2.type);

      bool passCuts = true;
      for(auto cut : stats.bset) {

        if(!passCuts) break;
        else if (cut == "DiscrByDeltaR") passCuts = passCuts && (part1.DeltaR(part2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(part1.Phi() - part2.Phi())), stats.pmap.at("CosDphiCut"));
        else if(cut == "DiscrByDeltaPt") passCuts = passCuts && passCutRange(part1.Pt() - part2.Pt(), stats.pmap.at("DeltaPtCutValue"));
        else if(cut == "DiscrByCDFzeta2D") {
          double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2).first
            + stats.dmap.at("PZetaVisCutCoefficient") * getPZeta(part1, part2).second;
          passCuts = passCuts && passCutRange(CDFzeta, stats.pmap.at("CDFzeta2DCutValue"));
        }
        else if(cut == "DiscrByDeltaPtDivSumPt") {
          double ptDiv = (part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt());
          passCuts = passCuts && passCutRange(ptDiv, stats.pmap.at("DeltaPtDivSumPtCutValue"));
        }
        else if (cut == "DiscrByMassReco") {
          double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
          passCuts = passCuts && passCutRange(diMass, stats.pmap.at("MassCut"));
        }
        else if(cut == "DiscrByCosDphiPtAndMet"){
          double CosDPhi1 = cos(absnormPhi(part1.Phi() - _MET->phi()));
          passCuts = passCuts && passCutRange(CosDPhi1, stats.pmap.at("CosDphiPtAndMetCut"));
        }
      	// ---------- New: DY Z' team ---------- //
      	else if(cut == "DiscrByCosDphiLeadPtAndMet"){
      	  if(part1.Pt() > part2.Pt()){
      	    llep = part1;
      	  }
      	  else{
      	    llep = part2;
      	  }
      	  double CosDPhiLead = cos(absnormPhi(llep.Phi() - _MET->phi()));
      	  passCuts = passCuts && passCutRange(CosDPhiLead, stats.pmap.at("CosDphiLeadPtAndMetCut"));
      	}
        else if (cut == "DiscrByAbsCosDphiLeadPtandMet"){
          if (part1.Pt() > part2.Pt()){
            llep = part1;
          }
          else{
            llep = part2;
          }
          double CosDPhiLead_forabs = cos(absnormPhi(llep.Phi() - _MET->phi()));
          passCuts = passCuts && passCutRangeAbs(CosDPhiLead_forabs, stats.pmap.at("AbsCosDphiLeadPtAndMetCut"));
        }
      	else if(cut == "DiscrByMtLeadPtAndMet"){
      	   if(part1.Pt() > part2.Pt()){
      	     llep = part1;
      	  }
      	  else{
      	    llep = part2;
      	  }
      	  double mTlead = calculateLeptonMetMt(llep);
      	  passCuts = passCuts && passCutRange(mTlead, stats.pmap.at("MtLeadPtAndMetCut"));
      	}
      	else if(cut == "DiscrByDiLepMassDeltaPt"){
      	  double dilepmass = CalculateDiLepMassDeltaPt(part1, part2, mass1, mass2);
      	  passCuts = passCuts && passCutRange(dilepmass, stats.pmap.at("DiLeadMassDeltaPtCut"));
      	}
      	// ---------------------------------------- //
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }  // end of for loop over bset.
      if (stats.bfind("DiscrByOSLSType")){
        //   if it is 1 or 0 it will end up in the bool std::map!!
        if(stats.bfind("DiscrByOSLSType") && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
      }else if (stats.dmap.find("DiscrByOSLSType") != stats.dmap.end() ){
        if(lep1.charge(i1) * lep2.charge(i2) > 0) continue;
      }else if (stats.smap.find("DiscrByOSLSType") != stats.smap.end() ){
        if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
        else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge(i1) * lep2.charge(i2) >= 0)) continue;
      }

      ///Particles that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2
      if(passCuts)
        active_part->at(ePosFin)->push_back(i1*BIG_NUM + i2);
    }
  }
}

// ---------- New function: DY Z' team ---------- //
double Analyzer::CalculateDiLepMassDeltaPt(const TLorentzVector& part1, const TLorentzVector& part2, const float mass1, const float mass2){

      double pt1  = part1.Pt();
      double eta1 = part1.Eta();
      double phi1 = part1.Phi();
      double pt2  = part2.Pt();
      double eta2 = part2.Eta();
      double phi2 = part2.Phi();
      double mass3= 0.0;

      double px1 = pt1*cos(phi1);
      double py1 = pt1*sin(phi1);
      double pz1 = pt1*sinh(eta1);
      double E1  = sqrt( pow(pt1*cosh(eta1),2) + pow(mass1,2) );
      double px2 = pt2*cos(phi2);
      double py2 = pt2*sin(phi2);
      double pz2 = pt2*sinh(eta2);
      double E2  = sqrt( pow(pt2*cosh(eta2),2) + pow(mass2,2) );
      double px3 = -(px1 + px2);
      double py3 = -(py1 + py2);
      double pz3 = -0.0;
      double E3  = sqrt( pow(px3,2) + pow(py3,2) + pow(pz3,2) + pow(mass3,2) );

      double diMass= sqrt(pow(E1+E2+E3,2)-pow(px1+px2+px3,2)-pow(py1+py2+py3,2)-pow(pz1+pz2+pz3,2));

      return diMass;
}
// ---------------------------------------- //

/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonJetCombos(Lepton& lep1, Jet& jet1, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);
  if(!lep1.needSyst(syst) && !jet1.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  TLorentzVector llep1, ljet1;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(ePos1)) {
    llep1 = lep1.p4(ij2);
    for(auto ij1 : *active_part->at(ePos2)) {
      ljet1 = _Jet->p4(ij1);

      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = (ljet1.DeltaR(llep1) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByDeltaEta") passCuts = passCutRange(abs(ljet1.Eta() - llep1.Eta()), stats.pmap.at("DeltaEtaCut"));
        else if(cut == "DiscrByDeltaPhi") passCuts = passCutRange(absnormPhi(ljet1.Phi() - llep1.Phi()), stats.pmap.at("DeltaPhiCut"));
        else if(cut == "DiscrByOSEta") passCuts = (ljet1.Eta() * llep1.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCutRange((ljet1+llep1).M(), stats.pmap.at("MassCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCutRange(cos(absnormPhi(ljet1.Phi() - llep1.Phi())), stats.pmap.at("CosDphiCut"));
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(ePosFin)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eDiJet)) return;
  std::string systname = syst_names.at(syst);
  if(systname!="orig"){
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eDiJet)=goodParts[CUTS::eDiJet];
      return;
    }
  }
  TLorentzVector jet1(0,0,0,0), jet2(0,0,0,0);
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(CUTS::eRJet2)) {
    jet2 = _Jet->p4(ij2);
    for(auto ij1 : *active_part->at(CUTS::eRJet1)) {
      if(ij1 == ij2) continue;
      jet1 = _Jet->p4(ij1);

      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = passCuts && (jet1.DeltaR(jet2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(jet1.Eta() - jet2.Eta()), stats.pmap.at("DeltaEtaCut"));
        else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(jet1.Phi() - jet2.Phi()), stats.pmap.at("DeltaPhiCut"));
        else if(cut == "DiscrByOSEta") passCuts = passCuts && (jet1.Eta() * jet2.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCuts && passCutRange((jet1+jet2).M(), stats.pmap.at("MassCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats.pmap.at("CosDphiCut"));
        else if(cut == "RejectForwardDiJetPairs") passCuts = passCuts && ( !passCutRangeAbs(jet1.Eta(), stats.pmap.at("ForwardEtaRange")) || !passCutRangeAbs(jet2.Eta(), stats.pmap.at("ForwardEtaRange")) );

        /*
        {
        	std::cout << "Checking if this is a forward dijet pair..." << std::endl;
        	std::cout << "Jet1 eta = " << jet1.Eta() << ", jet2 eta = " << jet2.Eta() << std::endl;
        	// Check the first jet of the pair:
        	std::cout << "Jet 1 in the forward eta range? " << passCutRangeAbs(jet1.Eta(), stats.pmap.at("ForwardEtaRange")) << std::endl;
        	// Then go onto the second jet of the pair:
        	passCuts = passCuts && !passCutRangeAbs(jet2.Eta(), stats.pmap.at("ForwardEtaRange"));
        	std::cout << "Jet 2 in the forward eta range? " << passCutRangeAbs(jet2.Eta(), stats.pmap.at("ForwardEtaRange")) << std::endl;
        	if(passCuts){ std::cout << "dijet pair not coming from both jets in the forward region" << std::endl;}
        	else{std::cout << "dijet pair -- coming -- from both jets in the forward region" << std::endl;}
        }
        */
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(CUTS::eDiJet)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}

///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(std::vector<TLorentzVector>::const_iterator lepit= lep.begin(); lepit != lep.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Calculates the Pzeta value
std::pair<double, double> Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + _MET->px();
  double py = visPy + _MET->py();
  return std::make_pair(px*zetaX + py*zetaY, visPx*zetaX + visPy*zetaY);
}

void Analyzer::checkParticleDecayList(){
  std::fstream file;
  file.open("particle_decay_list.txt", std::fstream::in | std::fstream::out);
  int s;
  if (file.is_open()){
    std::cout << "Warning, file already exists.";
    std::cout << "Do you wish to clear the file? 1 for yes; 0 for no.";
    std::cin >> s;
    if (s == 1)
      {file.open("particle_decay_list.txt", std::ios::out | std::ios::trunc);
	file.close();
	std::cout << "You have cleared the file.";}
  }
}

void Analyzer::writeParticleDecayList(int event){  //01.16.19
  BOOM->GetEntry(event);
  std::fstream file("particle_decay_list.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  if (file.is_open()){
  for (unsigned p=0; p < _Gen->size(); p++){
    file << std::setw(2) << p << std::setw(2) << " " << std::setw(8) << "pdg_id: " << std::setw(4) << abs(_Gen->pdg_id[p]) << std::setw(2) << "  " << std::setw(5) << "p_T: " << std::setw(10) << _Gen->pt(p) << std::setw(2) << "  " << std::setw(5) << "phi: " << std::setw(10) << _Gen->phi(p) << std::setw(10) << "status: " << std::setw(3) << _Gen->status[p] << std::setw(2) << "  " << std::setw(7) << "mind:" << std::setw(1) << " " << std::setw(2) << _Gen->genPartIdxMother[p] << "\n";
  }
  file << "----------" << "\n";
  file.close();}
  else std::cout << "Unable to open file." << std::endl;
  return;
}

std::multimap<int,int> Analyzer::readinJSON(std::string year){ //NEW:  function for reading in the JSON file in the format that I established.  This happens once before events start entering preprocess.
  // Newer: this function will take the year as an input to set the JSON file to be read.
  // The naming convention is json201X.txt, make sure your JSON files match this convention.
  std::fstream fs((PUSPACE+"json"+year+".txt").c_str(), std::fstream::in); // Take the info of the JSON file as input.
  std::string line; //NEW:  need this since we'll be analyzing line by line.
  std::multimap<int,int> json_line_dict; //NEW:  we'll work with the information as a multimap (C++ equivalent of a python dictionary).

  while(!fs.eof()){ //NEW:  while the file is open...
      getline(fs, line); //NEW:  grab each line.
      if(line.size() == 0) continue;  //NEW:  if there's nothing in the line, skip it.
      std::vector<std::string> vals = string_split(line, {","});  //NEW:  put the contents of the line into a vector.
      int run_num_line = (stringtotype<int>(vals[0])); //NEW:  the first element of the line is the run_number.
      for(size_t p=1; p < vals.size(); p=p+1){  //NEW:  step through the remaining elements of the line.
	json_line_dict.insert (std::make_pair(run_num_line,stringtotype<int>(vals[p]))); //NEW:  push the remaining lumisection bounds as each an element corresponding to the run_number key.
      }
  }
  fs.close(); // Close the file once the loop has gone over all lines in the file.

  return json_line_dict;  //NEW:  returns the multimap that you need for proper filtering.
}

double Analyzer::getZBoostWeight(){
  double boostweigth=1.;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1)){
    //cout<<" Z or W " <<std::endl;
    double boostz = 0;
    if(active_part->at(CUTS::eGZ)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
    }
    if(active_part->at(CUTS::eGW)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
    }
    if(boostz > 0 && boostz <= 50) {boostweigth = 1.1192;}
    else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1034;}
    else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0675;}
    else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0637;}
    else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0242;}
    else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9453;}
    else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8579;}
    else if (boostz >= 600) {boostweigth = 0.7822;}
    else {boostweigth = 1;}
  }
  return boostweigth;
}

double Analyzer::getZBoostWeightSyst(int ud){
  double boostweigth=1.;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->s\
ize() ==1)){
    double boostz = 0;
    if(active_part->at(CUTS::eGZ)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
    }
    if(active_part->at(CUTS::eGW)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
    }
    if (ud == 0){
      if(boostz > 0 && boostz <= 50) {boostweigth = 1.1192;}// 1.0942, 1.1192, 1.1442 5.26
      else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1034;}// 1.0901, 1.1034, 1.1167
      else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0675;}// 1.0559, 1.0675, 1.0791
      else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0637;}// 1.0511, 1.0637, 1.0763
      else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0242;}// 1.011, 1.0242, 1.0374
      else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9453;}// 0.9269, 0.9453, 0.9637
      else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8579;}// 0.8302, 0.8579, 0.8856
      else if (boostz >= 600) {boostweigth = 0.7822;}// 0.6692, 0.7822, 0.8952
      else {boostweigth = 1;}}

    else if (ud == -1){
      if(boostz > 0 && boostz <= 50) {boostweigth = 1.0942;}// 1.0942, 1.1192, 1.1442 5.26
      else if (boostz > 50 && boostz <= 100) {boostweigth = 1.0901;}// 1.0901, 1.1034, 1.1167
      else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0559;}// 1.0559, 1.0675, 1.0791
      else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0511;}// 1.0511, 1.0637, 1.0763
      else if (boostz > 200 && boostz <= 300) {boostweigth = 1.011;}// 1.011, 1.0242, 1.0374
      else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9269;}// 0.9269, 0.9453, 0.9637
      else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8302;}// 0.8302, 0.8579, 0.8856
      else if (boostz >= 600) {boostweigth = 0.6692;}// 0.6692, 0.7822, 0.8952
      else {boostweigth = 1;}}

    else if (ud == 1){
      if(boostz > 0 && boostz <= 50) {boostweigth = 1.1442;}// 1.0942, 1.1192, 1.1442 5.26
      else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1167;}// 1.0901, 1.1034, 1.1167
      else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0791;}// 1.0559, 1.0675, 1.0791
      else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0763;}// 1.0511, 1.0637, 1.0763
      else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0374;}// 1.011, 1.0242, 1.0374
      else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9637;}// 0.9269, 0.9453, 0.9637
      else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8856;}// 0.8302, 0.8579, 0.8856
      else if (boostz >= 600) {boostweigth = 0.8952;}// 0.6692, 0.7822, 0.8952
      else {boostweigth = 1;}}

  }
  return boostweigth;
}

double Analyzer::getWkfactor(){
  double kfactor=1.;
  if(!isWSample)
    return kfactor;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGW)->size() ==1)){
    //this k-factor is not computed for the low mass W!
    double wmass=_Gen->p4(active_part->at(CUTS::eGW)->at(0)).M();
    if(wmass<100){
      return 1.;
    }
    if(active_part->at(CUTS::eGTau)->size()){
      kfactor=k_ele_h->GetBinContent(k_ele_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGMuon)->size()){
      kfactor=k_mu_h->GetBinContent(k_mu_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGElec)->size()){
      kfactor=k_tau_h->GetBinContent(k_tau_h->FindBin(wmass));
    }
  }
  return kfactor;
}

// The function below applies the Z-pT corrections derived by the SUSY PAG - Ref. AN2015_267_v10.pdf
double Analyzer::getZpTWeight() {

  double zPtBoost = 1.;

  if(!((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1
    && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1))) return zPtBoost;


  double zMass = 0, zPT = 0;

  if(active_part->at(CUTS::eGZ)->size() == 1) {
      zMass = _Gen->mass(active_part->at(CUTS::eGZ)->at(0));
      zPT = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
  }
  if(active_part->at(CUTS::eGW)->size() == 1) {
    zMass = _Gen->mass(active_part->at(CUTS::eGW)->at(0));
    zPT = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
  }

  if(20 <= zMass && zMass < 60) {
     if(20 <= zPT && zPT < 40) zPtBoost = 1.04;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.12;
     else if(60 <= zPT && zPT < 80) zPtBoost = 1.09;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.09;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.14;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.29;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.21;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.47;
     else if(220 <= zPT && zPT < 300) zPtBoost = 0.78;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 1.56;
  }

  else if(60 <= zMass && zMass < 120) {
     if(20 <= zPT && zPT < 40) zPtBoost = 0.99;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.05;
     else if(60 <= zPT && zPT < 80) zPtBoost = 1.12;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.16;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.17;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.17;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.17;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.18;
     else if(220 <= zPT && zPT < 300) zPtBoost = 1.13;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 1.03;
  }

  else if(120 <= zMass && zMass < 160) {
     if(20 <= zPT && zPT < 40) zPtBoost = 1.07;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.13;
     else if(60 <= zPT && zPT < 80) zPtBoost = 1.16;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.21;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.22;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.27;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.28;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.17;
     else if(220 <= zPT && zPT < 300) zPtBoost = 1.35;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 1.06;
  }

  else if(160 <= zMass && zMass < 200) {
     if(20 <= zPT && zPT < 40) zPtBoost = 1.17;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.15;
     else if(60 <= zPT && zPT < 80) zPtBoost = 1.21;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.09;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.35;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.21;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.42;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.43;
     else if(220 <= zPT && zPT < 300) zPtBoost = 1.28;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 1.09;
  }

  else if(200 <= zMass && zMass < 240) {
     if(20 <= zPT && zPT < 40) zPtBoost = 1.22;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.21;
     else if(60 <= zPT && zPT < 80) zPtBoost = 0.97;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.42;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.41;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.17;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.30;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.21;
     else if(220 <= zPT && zPT < 300) zPtBoost = 1.45;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 0.85;
  }

  else if(240 <= zMass && zMass < 10000) {
     if(20 <= zPT && zPT < 40) zPtBoost = 1.24;
     else if(40 <= zPT && zPT < 60) zPtBoost = 1.47;
     else if(60 <= zPT && zPT < 80) zPtBoost = 1.26;
     else if(80 <= zPT && zPT < 100) zPtBoost = 1.48;
     else if(100 <= zPT && zPT < 120) zPtBoost = 1.61;
     else if(120 <= zPT && zPT < 140) zPtBoost = 1.28;
     else if(140 <= zPT && zPT < 180) zPtBoost = 1.44;
     else if(180 <= zPT && zPT < 220) zPtBoost = 1.41;
     else if(220 <= zPT && zPT < 300) zPtBoost = 1.47;
     else if(300 <= zPT && zPT < 10000) zPtBoost = 0.92;
  }

  //std::cout << "V mass = " << zMass << ", V pt = " << zPT << ", zPtBoost = " << zPtBoost << std::endl;

  return zPtBoost;
}

// These are weights derived by the VBF SUSY team - Run II (Kyungmin Park)
double Analyzer::getZpTWeight_vbfSusy(std::string year) {
    double zPtBoost = 1.;

    // std::cout << "Year (zboostwgt) = " << year << std::endl;
    if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1)){
      double zPT = 0;

      if(active_part->at(CUTS::eGZ)->size() ==1) {
          zPT = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
      }
      if(active_part->at(CUTS::eGW)->size() ==1) {
          zPT = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
      }

    // Correction factors derived from 0-lepton channel Zjets CR (2016)
    if(year.compare("2016") == 0){
      // std::cout << "This is 2016" << std::endl;
      if (0 <= zPT && zPT < 20) zPtBoost = 1.00;
      else if(20 <= zPT && zPT < 40) zPtBoost = 0.98;
      else if(40 <= zPT && zPT < 60) zPtBoost = 1.02;
      else if(60 <= zPT && zPT < 80) zPtBoost = 1.06;
      else if(80 <= zPT && zPT < 100) zPtBoost = 1.05;
      else if(100 <= zPT && zPT < 120) zPtBoost = 1.01;
      else if(120 <= zPT && zPT < 140) zPtBoost = 0.99;
      else if(140 <= zPT && zPT < 160) zPtBoost = 0.98;
      else if(160 <= zPT && zPT < 180) zPtBoost = 0.96;
      else if(180 <= zPT && zPT < 200) zPtBoost = 0.93;
      else if(200 <= zPT && zPT < 220) zPtBoost = 0.92;
      else if(220 <= zPT && zPT < 240) zPtBoost = 0.86;
      else if(240 <= zPT && zPT < 260) zPtBoost = 0.87;
      else if(260 <= zPT && zPT < 280) zPtBoost = 0.85;
      else if(280 <= zPT && zPT < 300) zPtBoost = 0.86;
      else if(300 <= zPT && zPT < 5000) zPtBoost = 0.88;
    }

    // Correction factors derived from 0-lepton channel Zjets CR (2017+2018)
    if(year.compare("2017") == 0 || year.compare("2018") == 0){
      // std::cout << "This is 2017 or 2018" << std::endl;
      if (0 <= zPT && zPT < 20) zPtBoost = 0.94;
      else if(20 <= zPT && zPT < 40) zPtBoost = 1.21;
      else if(40 <= zPT && zPT < 60) zPtBoost = 1.20;
      else if(60 <= zPT && zPT < 80) zPtBoost = 1.19;
      else if(80 <= zPT && zPT < 100) zPtBoost = 1.14;
      else if(100 <= zPT && zPT < 120) zPtBoost = 1.05;
      else if(120 <= zPT && zPT < 140) zPtBoost = 1.01;
      else if(140 <= zPT && zPT < 160) zPtBoost = 0.99;
      else if(160 <= zPT && zPT < 180) zPtBoost = 0.96;
      else if(180 <= zPT && zPT < 200) zPtBoost = 0.94;
      else if(200 <= zPT && zPT < 220) zPtBoost = 0.90;
      else if(220 <= zPT && zPT < 240) zPtBoost = 0.89;
      else if(240 <= zPT && zPT < 260) zPtBoost = 0.86;
      else if(260 <= zPT && zPT < 280) zPtBoost = 0.86;
      else if(280 <= zPT && zPT < 300) zPtBoost = 0.84;
      else if(300 <= zPT && zPT < 5000) zPtBoost = 0.83;
    }
  }

    return zPtBoost;
}

// Brenda FE -- additional PU weights that correct the NPV distribution
double Analyzer::getNPVWeight(bool applyWeights, float npv){

  double npv_weight = 1.0;
  if(!applyWeights) return npv_weight;

  npv_weight = histnpvwgt->GetBinContent(histnpvwgt->FindBin(npv));
  // std::cout << "NPV = " << npv << ", bin in NPV histo = " << (histnpvwgt->FindBin(npv)) << ", NPV weight = " << npv_weight << std::endl;
  return npv_weight;
}

////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram(std::string year) {

  if(!isData && distats["Run"].bfind("ApplyGenWeight") && gen_weight == 0.0) return;

  if(isData && blinded && maxCut == SignalRegion) return;

  const std::vector<std::string>* groups = histo.get_groups();

  if(!isData){
    wgt = 1.;
    //wgt *= getTopBoostWeight(); //01.15.19
    if(distats["Run"].bfind("UsePileUpWeight")) wgt*= pu_weight;
    if(distats["Run"].bfind("ApplyGenWeight")) wgt *= (gen_weight > 0) ? 1.0 : -1.0;
    //add weight here
    //if(distats["Run"].bfind("ApplyTauIDSF")) wgt *= getTauIdDataMCScaleFactor(0);
    // Arguments: getTauIdSFs(bool getTauIDsf, bool getTauIDbyDMsfs, bool getAntiElesf, bool getAntiMusf, std::string uncertainty)
    // std::cout << "ApplyTauIDSF = " << distats["Run"].bfind("ApplyTauIDSF") << ", TauIdSFsByDM = " << distats["Run"].bfind("TauIdSFsByDM") << ", ApplyTauAntiEleSF = " << distats["Run"].bfind("ApplyTauAntiEleSF") << ", ApplyTauAntiMuSF = " << distats["Run"].bfind("ApplyTauAntiMuSF") << std::endl;
    if(distats["Run"].bfind("ApplyTauIDSF")) wgt *= getTauIdSFs(true, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("ApplyTauAntiEleSF"), distats["Run"].bfind("ApplyTauAntiMuSF"), "");

    // Apply Z-boost weights derived for ISR+stau analysis (SUS-19-002)
    if(distats["Run"].bfind("ApplyISRZBoostSF") && isVSample){
      //wgt *= getZBoostWeight();
      wgt *= getZBoostWeightSyst(0);
      boosters[0] = getZBoostWeightSyst(0); //06.02.20
      boosters[1] = getZBoostWeightSyst(-1);  //06.02.20
      boosters[2] = getZBoostWeightSyst(1);  //06.02.20
    }
    else if(distats["Run"].bfind("ApplySUSYZBoostSF") && isVSample){
      wgt *= getZpTWeight();
    }
    else if(distats["Run"].bfind("ApplyVBFSusyZBoostSF") && isVSample){
      wgt *= getZpTWeight_vbfSusy(year);
    }
    if(distats["Run"].bfind("ApplyWKfactor")){
      wgt *= getWkfactor();
    }
    if(distats["Run"].bfind("ApplyNPVWeight")){
      wgt *= getNPVWeight(true, bestVertices);
    }

    // Apply L1 prefiring weights for 2016/2017
    if(distats["Run"].bfind("ApplyL1PrefiringWeight")){ // January 20, 2021 - Brenda FE
      // std::cout << "L1 prefiring weight = " << l1prefiringwgt << std::endl;
      wgt *= l1prefiringwgt; // nominal value

    }

    wgt *= getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]); //01.16.19
  }else  wgt=1.;
  //backup current weight
  backup_wgt=wgt;


  for(size_t i = 0; i < syst_names.size(); i++) {

    for(Particle* ipart: allParticles) ipart->setCurrentP(i);
    _MET->setCurrentP(i);

    active_part = &syst_parts.at(i);

    if(i == 0) {
      active_part = &goodParts;
      fillCuts(true);
      if(isSignalMC && !finalInputSignal) continue;
      for(auto it: *groups) {
        fill_Folder(it, maxCut, histo, false);
      }
      //if(!fillCuts(false)) {
      //  fill_Tree();
      //}

    }else{
      wgt=backup_wgt;

      if(syst_names[i].find("weight")!=std::string::npos){
        // ---------- Tau ID scale factors ---------- //
        if(syst_names[i]=="Tau_weight_Up"){
          if(distats["Run"].bfind("ApplyTauIDSF")) {
            // Arguments: getTauIdSFs(bool getTauIDsf, bool getTauIDbyDMsfs, bool getAntiElesf, bool getAntiMusf, std::string uncertainty)
            wgt /= getTauIdSFs(true, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("ApplyTauAntiEleSF"), distats["Run"].bfind("ApplyTauAntiMuSF"), "");
            wgt *= getTauIdSFs(true, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("ApplyTauAntiEleSF"), distats["Run"].bfind("ApplyTauAntiMuSF"), "Up");
          }
        }else if(syst_names[i]=="Tau_weight_Down"){
          if(distats["Run"].bfind("ApplyTauIDSF")) {
            wgt /= getTauIdSFs(true, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("ApplyTauAntiEleSF"), distats["Run"].bfind("ApplyTauAntiMuSF"), "");
            wgt *= getTauIdSFs(true, distats["Run"].bfind("TauIdSFsByDM"), distats["Run"].bfind("ApplyTauAntiEleSF"), distats["Run"].bfind("ApplyTauAntiMuSF"), "Down");
          }
        }
        // ---------- Pileup weights ---------- //
        if(syst_names[i]=="Pileup_weight_Up"){
          if(distats["Run"].bfind("UsePileUpWeight")) {
            wgt/=   pu_weight;
            // wgt *=  hPU_up[(int)(nTruePU+1)];
            wgt *= hist_pu_wgt_up->GetBinContent(hist_pu_wgt_up->FindBin(nTruePU));
          }
        }else if(syst_names[i]=="Pileup_weight_Down"){
          if(distats["Run"].bfind("UsePileUpWeight")) {
            wgt/=   pu_weight;
            // wgt *=  hPU_down[(int)(nTruePU+1)];
            wgt *= hist_pu_wgt_do->GetBinContent(hist_pu_wgt_do->FindBin(nTruePU));
          }
        }
        // --------- Prefiring weights ----------- //
        if(syst_names[i]=="L1Prefiring_weight_Up"){
          if(distats["Run"].bfind("ApplyL1PrefiringWeight")){
            wgt /= l1prefiringwgt;
            wgt *= l1prefiringwgt_up;
          }
        } else if(syst_names[i]=="L1Prefiring_weight_Down"){
          if(distats["Run"].bfind("ApplyL1PrefiringWeight")){
            wgt /= l1prefiringwgt;
            wgt *= l1prefiringwgt_dn;
          }
        }
      }

      if(syst_names[i].find("Btag")!=std::string::npos){ //01.16.19
        if(syst_names[i]=="Btag_Up"){
          wgt/=getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]);
          wgt*=getBJetSFResUp(CUTS::eRBJet, _Jet->pstats["BJet"]);
        }else if(syst_names[i]=="Btag_Down"){
          wgt/=getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]);
          wgt*=getBJetSFResDown(CUTS::eRBJet, _Jet->pstats["BJet"]);
        }
      }

      ///---06.06.20---
      if(syst_names[i].find("ISR_weight")!=std::string::npos){ //07.09.18
        if(syst_names[i]=="ISR_weight_up"){
          if(distats["Run"].bfind("ApplyISRZBoostSF") && isVSample) {
            wgt/=boosters[0];
            wgt*=boosters[2];
          }
        }else if(syst_names[i]=="ISR_weight_down"){
          if(distats["Run"].bfind("ApplyISRZBoostSF") && isVSample) {
            wgt/=boosters[0];
            wgt*=boosters[1];
          }
        }
      }
      ///---06.02.20
    }

    //get the non particle conditions:
    for(auto itCut : nonParticleCuts){
      active_part->at(itCut)=goodParts.at(itCut);
    }

    if(!fillCuts(false)) continue;
    if(isSignalMC && !finalInputSignal) continue; // Brenda FE, Apr 20, 2021

    for(auto it: *syst_histo.get_groups()) {
      fill_Folder(it, i, syst_histo, true);
    }
    wgt=backup_wgt;
  }

  for(Particle* ipart: allParticles) ipart->setCurrentP(0);
   _MET->setCurrentP(0);
  active_part = &goodParts;

}

///Function that fills up the histograms
void Analyzer::fill_Folder(std::string group, const int max, Histogramer &ihisto, bool issyst) {
  /* be aware in this function
   * the following definition is used:
   * histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
   * so each histogram knows the group, max and weight!
   */

  if(group == "FillRun" && (&ihisto==&histo)) {
    if(crbins != 1) {
      for(int i = 0; i < crbins; i++) {
        std::cout << "nominal, crbins !=1" << std::endl;
        ihisto.addVal(false, group, i, "Events", 1);
        if(distats["Run"].bfind("ApplyGenWeight")) {
          //put the weighted events in bin 3
          ihisto.addVal(2, group,i, "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
        ihisto.addVal(wgt, group, i, "Weight", 1);
      }
    }
    else{

      if(isSignalMC){ // For signal samples (unskimmed)
        if(finalInputSignal){

	    	  // Bin 1 will contain only the events that correspond to a certain signal mass point
          ihisto.addVal(false, group,ihisto.get_maxfolder(), "Events", 1);

          if(distats["Run"].bfind("ApplyGenWeight")) {

            if(finalInputSignal){
              // Bin 3 contains the number of raw signal MC events passing a particular cut, same as bin1 but without any weights.
              ihisto.addVal(2, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
            }
            // Bin 4 contains number of raw MC events for a given signal point for which the gen_weight is positive.
            ihisto.addVal(3, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
          }
        }
        // Bin 5 contains the full number of events analyzed in the signal sample (including all signal points)
        ihisto.addVal(4,group,ihisto.get_maxfolder(),"Events",1.0);
	    }
	    else if(!isSignalMC || isData ){ // For backgrounds or data.
        // Bin 1 contains the number of events analyzed in a given sample.
	    	ihisto.addVal(false, group,ihisto.get_maxfolder(), "Events", 1);

        if(distats["Run"].bfind("ApplyGenWeight")) {
          // Bin 3 contains the raw (unweighted) number of events (data or MC) passing a certain selection.
          ihisto.addVal(2, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
	    }
      ihisto.addVal(wgt, group, ihisto.get_maxfolder(), "Weight", 1);
      ihisto.addVal(nTruePU, group, ihisto.get_maxfolder(), "PUWeight", 1);
    }

    // In all cases: Bin 2 contains the weighted number of events passing a certain selection.
    histAddVal(true, "Events");
    histAddVal(totalVertices, "NTotalVertices");
    histAddVal(bestVertices, "NVertices");

  } else if(group == "FillRun" && issyst) {

    if(isSignalMC){ // For signal samples (unskimmed)
      if(finalInputSignal){

        // Bin 1 will contain only the events that correspond to a certain signal mass point
        syst_histo.addVal(false, group,max, "Events", 1);

        if(distats["Run"].bfind("ApplyGenWeight")) {

          if(finalInputSignal){
            // Bin 3 contains the number of raw signal MC events passing a particular cut, same as bin1 but without any weights.
            syst_histo.addVal(2, group,max, "Events", (gen_weight > 0) ? 1.0 : -1.0);
          }
          // Bin 4 contains number of raw MC events for a given signal point for which the gen_weight is positive.
          syst_histo.addVal(3, group,max, "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
      }
      // Bin 5 contains the full number of events analyzed in the signal sample (including all signal points)
      syst_histo.addVal(4,group,max,"Events",1.0);
    }
    else if(!isSignalMC){ // For backgrounds
      // Bin 1 contains the number of events analyzed in a given sample.
      syst_histo.addVal(false, group,max, "Events", 1);

      if(distats["Run"].bfind("ApplyGenWeight")) {
        // Bin 3 contains the raw (unweighted) number of events (data or MC) passing a certain selection.
        syst_histo.addVal(2, group,max, "Events", (gen_weight > 0) ? 1.0 : -1.0);
      }
    }
    syst_histo.addVal(wgt, group, max, "Weight", 1);
    syst_histo.addVal(nTruePU, group, max, "PUWeight", 1);

    // In all cases: Bin 2 contains the weighted number of events passing a certain selection.
    histAddVal(true, "Events");
    histAddVal(totalVertices, "NTotalVertices");
    histAddVal(bestVertices, "NVertices");

  } else if(!isData && group == "FillGen") {

    // Add these to the default histograms to perform quick checks if necessary.
    histAddVal(nTruePU, "PUNTrueInt");
    histAddVal(generatorht, "HT");
    histAddVal(gen_weight, "Weight");
    histAddVal(l1prefiringwgt, "L1PrefiringWeight");

    std::vector<int> goodGenLeptons;

    int nhadtau = 0;
    TLorentzVector genVec(0,0,0,0);
    int i = 0;
    for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it!=active_part->at(CUTS::eGTau)->end(); it++, i++) {
	//std::cout << "CUTS::eGTau #" << i << ", index in Gen = " << (*it) << ", status = " << _Gen->status[(*it)] << std::endl;
      goodGenLeptons.push_back(*it);

      int nu = active_part->at(CUTS::eGNuTau)->at(i);
      if(nu != -1) {
        genVec = _Gen->p4(*it) - _Gen->p4(nu);
        histAddVal(genVec.Pt(), "HadTauPt");
        histAddVal(genVec.Eta(), "HadTauEta");
        nhadtau++;
      }
      histAddVal(_Gen->energy(*it), "TauEnergy");
      histAddVal(_Gen->pt(*it), "TauPt");
      histAddVal(_Gen->eta(*it), "TauEta");
      histAddVal(_Gen->phi(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=active_part->at(CUTS::eGTau)->end(); it2++) {
        histAddVal(diParticleMass(_Gen->p4(*it),_Gen->p4(*it2), "none"), "DiTauMass");
      }
    }

    for(vec_iter genhadtau_it = active_part->at(CUTS::eGHadTau)->begin(); genhadtau_it != active_part->at(CUTS::eGHadTau)->end(); genhadtau_it++){
        histAddVal(_GenHadTau->pt(*genhadtau_it), "VisHadTauPt");
        histAddVal(_GenHadTau->eta(*genhadtau_it), "VisHadTauEta");
        histAddVal(_GenHadTau->phi(*genhadtau_it), "VisHadTauPhi");
        histAddVal(_GenHadTau->p4(*genhadtau_it).M(), "VisHadTauMass");
        histAddVal(_GenHadTau->decayMode[*genhadtau_it], "VisHadTauDecayMode");
    }
    histAddVal(active_part->at(CUTS::eGHadTau)->size(), "NVisHadTau");

    for(vec_iter matchedgentauh_it = active_part->at(CUTS::eGMatchedHadTau)->begin(); matchedgentauh_it != active_part->at(CUTS::eGMatchedHadTau)->end(); matchedgentauh_it++){
        histAddVal(_GenHadTau->pt(*matchedgentauh_it), "MatchedVisTauHPt");
        histAddVal(_GenHadTau->eta(*matchedgentauh_it), "MatchedVisTauHEta");
        histAddVal(_GenHadTau->phi(*matchedgentauh_it), "MatchedVisTauHPhi");
        histAddVal(_GenHadTau->p4(*matchedgentauh_it).M(), "MatchedVisTauHMass");
        histAddVal(_GenHadTau->decayMode[*matchedgentauh_it], "MatchedVisTauHDecayMode");
    }
    histAddVal(active_part->at(CUTS::eGHadTau)->size(), "NMatchedVisTauH");

    for(auto it : *active_part->at(CUTS::eGBJet)){
      histAddVal(_GenJet->pt(it), "BJPt");
      histAddVal(_GenJet->eta(it), "BJEta");
    }

    histAddVal(active_part->at(CUTS::eGTau)->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");
    histAddVal(active_part->at(CUTS::eGBJet)->size(), "NBJ");

    for(auto it : *active_part->at(CUTS::eGZ)) {
      histAddVal(_Gen->pt(it), "ZPt");
      histAddVal(_Gen->eta(it), "ZEta");
      histAddVal(_Gen->p4(it).M(), "ZMass");
    }
    histAddVal(active_part->at(CUTS::eGZ)->size(), "NZ");

    for(auto it : *active_part->at(CUTS::eGW)) {
      histAddVal(_Gen->pt(it), "WPt");
      histAddVal(_Gen->eta(it), "WEta");
      histAddVal(_Gen->p4(it).M(), "WMass");
    }
    histAddVal(active_part->at(CUTS::eGW)->size(), "NW");

    for(auto it : *active_part->at(CUTS::eGMuon)) {
      goodGenLeptons.push_back(it);
      histAddVal(_Gen->energy(it), "MuonEnergy");
      histAddVal(_Gen->pt(it), "MuonPt");
      histAddVal(_Gen->eta(it), "MuonEta");
      histAddVal(_Gen->phi(it), "MuonPhi");
    }
    histAddVal(active_part->at(CUTS::eGMuon)->size(), "NMuon");

    for(auto it: *active_part->at(CUTS::eGElec)){
      goodGenLeptons.push_back(it);
      histAddVal(_Gen->energy(it), "ElectronEnergy");
      histAddVal(_Gen->pt(it), "ElectronPt");
      histAddVal(_Gen->eta(it), "ElectronEta");
      histAddVal(_Gen->phi(it), "ElectronPhi");
    }
    histAddVal(active_part->at(CUTS::eGElec)->size(), "NElectron");

    histAddVal(gendilepmass, "ZDiLepMass");

    double mass=0;
    // int nb_leptons = 0;
    TLorentzVector lep1(0,0,0,0), lep2(0,0,0,0);
    /*
    for(size_t igen=0; igen<_Gen->size(); igen++){
      if(abs(_Gen->pdg_id[igen])==11 or abs(_Gen->pdg_id[igen])==13 or abs(_Gen->pdg_id[igen])==15){

	if( (abs(_Gen->pdg_id[igen])==11 or abs(_Gen->pdg_id[igen])==13) && _Gen->status[igen] != 1) continue;
	else if( (abs(_Gen->pdg_id[igen])==15) && _Gen->status[igen] != 2 ) continue;

	// std::cout << "Gen pdg_id = " << (abs(_Gen->pdg_id[igen])) << ", status = " << _Gen->status[igen] << std::endl;
    */
    for(size_t i=0; i<goodGenLeptons.size(); i++){
        int igen = goodGenLeptons.at(i);
	// std::cout << "Gen pdg_id = " << (abs(_Gen->pdg_id[igen])) << ", status = " << _Gen->status[igen] << std::endl;
       
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(igen);
          mass=(lep1+lep2).M();
          break;
        }else{
          lep1= _Gen->RecoP4(igen);
        }

        histAddVal(_Gen->p4(igen).Pt(), "LeptonPt");
        histAddVal(_Gen->p4(igen).Eta(), "LeptonEta");
        histAddVal(_Gen->p4(igen).Phi(), "LeptonPhi");
        histAddVal(_Gen->p4(igen).E(), "LeptonE");

        //nb_leptons++;
    }

    histAddVal(mass, "DiLeptonMass");
    histAddVal(goodGenLeptons.size(), "NLepton");

    for(size_t i_genjet=0; i_genjet < _GenJet->size(); i_genjet++){
      histAddVal(_GenJet->pt(i_genjet), "JetPt");
      histAddVal(_GenJet->eta(i_genjet), "JetEta");
      histAddVal(_GenJet->phi(i_genjet), "JetPhi");
      histAddVal(_GenJet->energy(i_genjet), "JetEnergy");
      histAddVal(_GenJet->genHadronFlavor[i_genjet], "JetHadronFlavor");
      histAddVal(_GenJet->genPartonFlavor[i_genjet], "JetPartonFlavor");
    }

    histAddVal(_GenJet->size(), "NJet");

    histAddVal(genmet_pt, "MetPt");
    histAddVal(genmet_phi, "MetPhi");

  } else if(group == "FillSpecial"){

    if(_Jet->size() > 0){
      if(index_minjmetdphi_formet > -1){
        histAddVal(_Jet->pt(index_minjmetdphi_formet), "InitialMinAbsDPhiMetJetPt");
        histAddVal(_Jet->eta(index_minjmetdphi_formet), "InitialMinAbsDPhiMetJetEta");
        histAddVal2(_Jet->eta(index_minjmetdphi_formet), _MET->pt(), "InitialMinAbsDPhiMetVsJetEta");
        histAddVal2(_Jet->eta(index_minjmetdphi_formet), _Jet->pt(index_minjmetdphi_formet), "InitialMinAbsDPhiMetJetPtvsEta");
        histAddVal(_Jet->phi(index_minjmetdphi_formet), "InitialMinAbsDPhiMetJetPhi");
        histAddVal(minDeltaPhiMet_formet, "InitialMinAbsDPhiMetJet");
        histAddVal2(minDeltaPhiMet_formet, _MET->pt(), "InitialMetvsMinAbsDPhiMetJet");
      }

      if(index_maxjmetdphi_formet > -1){
        histAddVal(_Jet->pt(index_maxjmetdphi_formet), "InitialMaxAbsDPhiMetJetPt");
        histAddVal(_Jet->eta(index_maxjmetdphi_formet), "InitialMaxAbsDPhiMetJetEta");
        histAddVal2(_Jet->eta(index_maxjmetdphi_formet), _MET->pt(), "InitialMaxAbsDPhiMetVsJetEta");
        histAddVal2(_Jet->eta(index_minjmetdphi_formet), _Jet->pt(index_maxjmetdphi_formet), "InitialMaxAbsDPhiMetJetPtvsEta");
        histAddVal(_Jet->phi(index_maxjmetdphi_formet), "InitialMaxAbsDPhiMetJetPhi");
        histAddVal(maxDeltaPhiMet_formet, "InitialMaxAbsDPhiMetJet");
        histAddVal2(maxDeltaPhiMet_formet, _MET->pt(), "InitialMetvsMaxAbsDPhiMetJet");
      }

      if(index_maxjetptprojonmet_minus_formet > -1){
        histAddVal(_Jet->pt(index_maxjetptprojonmet_minus_formet), "InitialMaxJetProjMetMinusPt");
        histAddVal(_Jet->eta(index_maxjetptprojonmet_minus_formet), "InitialMaxJetProjMetMinusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_minus_formet), _MET->pt(), "InitialMaxJetProjMetMinusVsEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_minus_formet), _Jet->pt(index_maxjetptprojonmet_minus_formet), "InitialMaxJetProjMetMinusPtvsEta");
        histAddVal(_Jet->phi(index_maxjetptprojonmet_minus_formet), "InitialMaxJetProjMetMinusPhi");
        histAddVal(maxjetptprojonmet_minus_formet, "InitialMaxJetProjMetMinusMag");
        histAddVal2(maxjetptprojonmet_minus_formet, _MET->pt(), "InitialMetvsMaxJetProjMetMinusMag");
      }

      if(index_maxjetptprojonmet_plus_formet > -1){
        histAddVal(_Jet->pt(index_maxjetptprojonmet_plus_formet), "InitialMaxJetProjMetPlusPt");
        histAddVal(_Jet->eta(index_maxjetptprojonmet_plus_formet), "InitialMaxJetProjMetPlusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_plus_formet), _MET->pt(), "InitialMaxJetProjMetPlusVsEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_plus_formet), _Jet->pt(index_maxjetptprojonmet_plus_formet), "InitialMaxJetProjMetPlusPtvsEta");
        histAddVal(_Jet->phi(index_maxjetptprojonmet_plus_formet), "InitialMaxJetProjMetPlusPhi");
        histAddVal(maxjetptprojonmet_plus_formet, "InitialMaxJetProjMetPlusMag");
        histAddVal2(maxjetptprojonmet_plus_formet, _MET->pt(), "InitialMetvsMaxJetProjMetPlusMag");
      }

      if(index_maxjetptprojonmet_plus_formet > -1 || index_maxjetptprojonmet_minus_formet > -1){
        int maxprojonmet_idx_formet = maxjetptprojonmet_minus_formet > maxjetptprojonmet_plus_formet ? index_maxjetptprojonmet_minus_formet : index_maxjetptprojonmet_plus_formet;
        float maxprojonmet_formet = std::max(maxjetptprojonmet_minus_formet, maxjetptprojonmet_plus_formet);

        histAddVal(_Jet->pt(maxprojonmet_idx_formet), "InitialMaxJetProjMetPt");
        histAddVal(_Jet->eta(maxprojonmet_idx_formet), "InitialMaxJetProjMetEta");
        histAddVal2(_Jet->eta(maxprojonmet_idx_formet), _MET->pt(), "InitialMaxJetProjMetVsEta");
        histAddVal2(_Jet->eta(maxprojonmet_idx_formet), _Jet->pt(maxprojonmet_idx_formet), "InitialMaxJetProjMetPtvsEta");
        histAddVal(_Jet->phi(maxprojonmet_idx_formet), "InitialMaxJetProjMetPhi");
        histAddVal(maxprojonmet_formet, "InitialMaxJetProjMetMag");
        histAddVal2(maxprojonmet_formet, _MET->pt(), "InitialMetvsMaxJetProjMetMag");
        histAddVal(normPhi(_Jet->phi(maxprojonmet_idx_formet) - _MET->phi()), "InitialMaxJetProjMetDeltaPhi");
        histAddVal2(normPhi(_Jet->phi(maxprojonmet_idx_formet) - _MET->phi()), _MET->pt(), "InitialMetvsMaxJetProjMetDeltaPhi");
        if( cos( normPhi(_Jet->phi(maxprojonmet_idx_formet) - _MET->phi()) ) > 0 ){
          histAddVal(1, "InitialMaxJetProjMetSign");
          histAddVal2(1, _MET->pt(), "InitialMetvsMaxJetProjMetSign");
        } else {
          histAddVal(-1, "InitialMaxJetProjMetSign");
          histAddVal2(-1, _MET->pt(), "InitialMetvsMaxJetProjMetSign");
        }
        histAddVal(cos( normPhi(_Jet->phi(maxprojonmet_idx_formet) - _MET->phi()) ), "InitialMaxJetProjMetCosDPhi");
      }
    }

    // Plots for projection of bad jets in EE noise region for 2017, loop over all jets
    if(additionalEEnoiseJets.size() > 0 ){

      float maxjetptprojonmet = 0.0, minjetptprojonmet = 99999.9;
      int index_maxjetptprojonmet = -1, index_minjetptprojonmet = -1;
      float maxjetptprojperpmet = 0.0, minjetptprojperpmet = 99999.9;
      int index_maxjetptprojperpmet = -1, index_minjetptprojperpmet = -1;


      for(size_t i=0; i < additionalEEnoiseJets.size(); i++){

        // std::cout << "Jet #" << additionalEEnoiseJets.at(i) << ", pt = " << _Jet->pt(additionalEEnoiseJets.at(i)) << ", eta = " << _Jet->eta(additionalEEnoiseJets.at(i)) << std::endl;
        TLorentzVector jetP4 = _Jet->p4(additionalEEnoiseJets.at(i));

        float jetrawpt = _Jet->RecoP4(additionalEEnoiseJets.at(i)).Pt() * ( 1.0 - _Jet->rawFactor[additionalEEnoiseJets.at(i)] );
        float deltaPhiMet = normPhi(jetP4.Phi() - _MET->phi());
        float projonmetaxis = jetP4.Pt() * cos(deltaPhiMet);
        float projperptometaxis = jetP4.Pt() * sin(deltaPhiMet);
        float ratioprojonmet = projonmetaxis / _MET->pt();
        float ratioprojperpmet = projperptometaxis / _MET->pt();
        // std::cout << "Jet #" << additionalEEnoiseJets.at(i) << ", pt = " << _Jet->pt(additionalEEnoiseJets.at(i)) << ", eta = " << _Jet->eta(additionalEEnoiseJets.at(i)) << ", raw pt = " << jetrawpt <<  std::endl;

        if(abs(projonmetaxis) > maxjetptprojonmet){
          maxjetptprojonmet = abs(projonmetaxis);
          index_maxjetptprojonmet = additionalEEnoiseJets.at(i);
        }

        if(abs(projonmetaxis) < minjetptprojonmet){
          minjetptprojonmet = abs(projonmetaxis);
          index_minjetptprojonmet = additionalEEnoiseJets.at(i);
        }

        if(abs(projperptometaxis) > maxjetptprojperpmet){
          maxjetptprojperpmet = abs(projperptometaxis);
          index_maxjetptprojperpmet = additionalEEnoiseJets.at(i);
        }

        if(abs(projperptometaxis) < minjetptprojperpmet){
          minjetptprojperpmet = abs(projperptometaxis);
          index_minjetptprojperpmet = additionalEEnoiseJets.at(i);
        }

        histAddVal(jetP4.Pt(), "EENoiseBadJetPt");
        histAddVal(jetrawpt, "EENoiseBadJetRawPt");
        histAddVal(jetP4.Eta(), "EENoiseBadJetEta");
        histAddVal(jetP4.Phi(), "EENoiseBadJetPhi");
        histAddVal(deltaPhiMet, "EENoiseBadJetMetDPhi");
        histAddVal(abs(deltaPhiMet), "EENoiseBadJetMetAbsDPhi");

        if(projonmetaxis < 0.0){
          histAddVal(-1, "EENoiseBadJetPtProjonMetSign");
        } else{
          histAddVal(1, "EENoiseBadJetPtProjonMetSign");
        }

        if(projperptometaxis < 0.0){
          histAddVal(-1, "EENoiseBadJetPtProjPerpMetSign");
        } else{
          histAddVal(1, "EENoiseBadJetPtProjPerpMetSign");
        }

        histAddVal(cos(deltaPhiMet), "EENoiseBadJetMetCosDPhi");
        histAddVal(projonmetaxis, "EENoiseBadJetPtProjonMet");
        histAddVal(abs(projonmetaxis), "EENoiseBadJetPtProjonMetMag");
        histAddVal(ratioprojonmet, "EENoiseBadJetPtProjonMetRatio");
        histAddVal(abs(ratioprojonmet), "EENoiseBadJetPtProjonMetRatioMag");

        histAddVal(sin(deltaPhiMet), "EENoiseBadJetMetSinDPhi");
        histAddVal(projperptometaxis, "EENoiseBadJetPtProjPerpMet");
        histAddVal(abs(projperptometaxis), "EENoiseBadJetPtProjPerpMetMag");
        histAddVal(ratioprojperpmet, "EENoiseBadJetPtProjPerpMetRatio");
        histAddVal(abs(ratioprojperpmet), "EENoiseBadJetPtProjPerpMetRatioMag");

      }

      histAddVal(additionalEEnoiseJets.size(), "EENoiseBadNJets");

      if(index_maxjetptprojonmet != -1){

        TLorentzVector maxprojonmetjetp4 = _Jet->p4(index_maxjetptprojonmet);

        float jetrawpt = _Jet->RecoP4(index_maxjetptprojonmet).Pt() * ( 1.0 - _Jet->rawFactor[index_maxjetptprojonmet] );
        float deltaPhiMet = normPhi(maxprojonmetjetp4.Phi() - _MET->phi());
        float projonmetaxis = maxprojonmetjetp4.Pt() * cos(deltaPhiMet);
        float ratioprojonmet = projonmetaxis / _MET->pt();

        histAddVal(maxprojonmetjetp4.Pt(), "MaxProjonMetEENoiseBadJetPt");
        histAddVal(jetrawpt, "MaxProjonMetEENoiseBadJetRawPt");
        histAddVal(maxprojonmetjetp4.Eta(), "MaxProjonMetEENoiseBadJetEta");
        histAddVal(maxprojonmetjetp4.Phi(), "MaxProjonMetEENoiseBadJetPhi");
        histAddVal(deltaPhiMet, "MaxProjonMetEENoiseBadJetMetDPhi");
        histAddVal(abs(deltaPhiMet), "MaxProjonMetEENoiseBadJetMetAbsDPhi");

        if(projonmetaxis < 0.0){
          histAddVal(-1, "MaxProjonMetEENoiseBadJetPtProjonMetSign");
        } else{
          histAddVal(1, "MaxProjonMetEENoiseBadJetPtProjonMetSign");
        }

        histAddVal(cos(deltaPhiMet), "MaxProjonMetEENoiseBadJetMetCosDPhi");
        histAddVal(projonmetaxis, "MaxProjonMetEENoiseBadJetPtProjonMet");
        histAddVal(abs(projonmetaxis), "MaxProjonMetEENoiseBadJetPtProjonMetMag");
        histAddVal(ratioprojonmet, "MaxProjonMetEENoiseBadJetPtProjonMetRatio");
        histAddVal(abs(ratioprojonmet), "MaxProjonMetEENoiseBadJetPtProjonMetRatioMag")

      }

      if(index_minjetptprojonmet != -1){
        TLorentzVector minprojonmetjetp4 = _Jet->p4(index_minjetptprojonmet);

        float jetrawpt = _Jet->RecoP4(index_minjetptprojonmet).Pt() * ( 1.0 - _Jet->rawFactor[index_minjetptprojonmet] );
        float deltaPhiMet = normPhi(minprojonmetjetp4.Phi() - _MET->phi());
        float projonmetaxis = minprojonmetjetp4.Pt() * cos(deltaPhiMet);
        float ratioprojonmet = projonmetaxis / _MET->pt();

        histAddVal(minprojonmetjetp4.Pt(), "MinProjonMetEENoiseBadJetPt");
        histAddVal(jetrawpt, "MinProjonMetEENoiseBadJetRawPt");
        histAddVal(minprojonmetjetp4.Eta(), "MinProjonMetEENoiseBadJetEta");
        histAddVal(minprojonmetjetp4.Phi(), "MinProjonMetEENoiseBadJetPhi");
        histAddVal(deltaPhiMet, "MinProjonMetEENoiseBadJetMetDPhi");
        histAddVal(abs(deltaPhiMet), "MinProjonMetEENoiseBadJetMetAbsDPhi");

        if(projonmetaxis < 0.0){
          histAddVal(-1, "MinProjonMetEENoiseBadJetPtProjonMetSign");
        } else{
          histAddVal(1, "MinProjonMetEENoiseBadJetPtProjonMetSign");
        }

        histAddVal(cos(deltaPhiMet), "MinProjonMetEENoiseBadJetMetCosDPhi");
        histAddVal(projonmetaxis, "MinProjonMetEENoiseBadJetPtProjonMet");
        histAddVal(abs(projonmetaxis), "MinProjonMetEENoiseBadJetPtProjonMetMag");
        histAddVal(ratioprojonmet, "MinProjonMetEENoiseBadJetPtProjonMetRatio");
        histAddVal(abs(ratioprojonmet), "MinProjonMetEENoiseBadJetPtProjonMetRatioMag")
      }

      if(index_maxjetptprojperpmet != -1){
        TLorentzVector maxprojperpmetjetp4 = _Jet->p4(index_maxjetptprojperpmet);

        float jetrawpt = _Jet->RecoP4(index_maxjetptprojperpmet).Pt() * ( 1.0 - _Jet->rawFactor[index_maxjetptprojperpmet] );
        float deltaPhiMet = normPhi(maxprojperpmetjetp4.Phi() - _MET->phi());
        float projperpmetaxis = maxprojperpmetjetp4.Pt() * sin(deltaPhiMet);
        float ratioprojperpmet = projperpmetaxis / _MET->pt();

        histAddVal(maxprojperpmetjetp4.Pt(), "MaxProjPerpMetEENoiseBadJetPt");
        histAddVal(jetrawpt, "MaxProjPerpMetEENoiseBadJetRawPt");
        histAddVal(maxprojperpmetjetp4.Eta(), "MaxProjPerpMetEENoiseBadJetEta");
        histAddVal(maxprojperpmetjetp4.Phi(), "MaxProjPerpMetEENoiseBadJetPhi");
        histAddVal(deltaPhiMet, "MaxProjPerpMetEENoiseBadJetMetDPhi");
        histAddVal(abs(deltaPhiMet), "MaxProjPerpMetEENoiseBadJetMetAbsDPhi");

        if(projperpmetaxis < 0.0){
          histAddVal(-1, "MaxProjPerpMetEENoiseBadJetPtProjonMetSign");
        } else{
          histAddVal(1, "MaxProjPerpMetEENoiseBadJetPtProjonMetSign");
        }

        histAddVal(sin(deltaPhiMet), "MaxProjPerpMetEENoiseBadJetMetSinDPhi");
        histAddVal(projperpmetaxis, "MaxProjPerpMetEENoiseBadJetPtProjonMet");
        histAddVal(abs(projperpmetaxis), "MaxProjPerpMetEENoiseBadJetPtProjonMetMag");
        histAddVal(ratioprojperpmet, "MaxProjPerpMetEENoiseBadJetPtProjonMetRatio");
        histAddVal(abs(ratioprojperpmet), "MaxProjPerpMetEENoiseBadJetPtProjonMetRatioMag")

      }

      if(index_minjetptprojperpmet != -1){
        TLorentzVector minprojperpmetjetp4 = _Jet->p4(index_minjetptprojperpmet);

        float jetrawpt = _Jet->RecoP4(index_minjetptprojperpmet).Pt() * ( 1.0 - _Jet->rawFactor[index_minjetptprojperpmet] );
        float deltaPhiMet = normPhi(minprojperpmetjetp4.Phi() - _MET->phi());
        float projperpmetaxis = minprojperpmetjetp4.Pt() * sin(deltaPhiMet);
        float ratioprojperpmet = projperpmetaxis / _MET->pt();

        histAddVal(minprojperpmetjetp4.Pt(), "MinProjPerpMetEENoiseBadJetPt");
        histAddVal(jetrawpt, "MinProjPerpMetEENoiseBadJetRawPt");
        histAddVal(minprojperpmetjetp4.Eta(), "MinProjPerpMetEENoiseBadJetEta");
        histAddVal(minprojperpmetjetp4.Phi(), "MinProjPerpMetEENoiseBadJetPhi");
        histAddVal(deltaPhiMet, "MinProjPerpMetEENoiseBadJetMetDPhi");
        histAddVal(abs(deltaPhiMet), "MinProjPerpMetEENoiseBadJetMetAbsDPhi");

        if(projperpmetaxis < 0.0){
          histAddVal(-1, "MinProjPerpMetEENoiseBadJetPtProjonMetSign");
        } else{
          histAddVal(1, "MinProjPerpMetEENoiseBadJetPtProjonMetSign");
        }

        histAddVal(sin(deltaPhiMet), "MinProjPerpMetEENoiseBadJetMetSinDPhi");
        histAddVal(projperpmetaxis, "MinProjPerpMetEENoiseBadJetPtProjonMet");
        histAddVal(abs(projperpmetaxis), "MinProjPerpMetEENoiseBadJetPtProjonMetMag");
        histAddVal(ratioprojperpmet, "MinProjPerpMetEENoiseBadJetPtProjonMetRatio");
        histAddVal(abs(ratioprojperpmet), "MinProjPerpMetEENoiseBadJetPtProjonMetRatioMag")
      }


    }


  } else if(fillInfo[group]->type == FILLER::Single) {
    Particle* part = fillInfo[group]->part;
    CUTS ePos = fillInfo[group]->ePos;

    // Added variable to calculate minimum deltaPhi between particle and MET
    float minDeltaPhiMet = 9999.9;
    float maxDeltaPhiMet = 0.0;
    float maxjetptprojonmet_plus = 0.0, maxjetptprojonmet_minus = 0.0;
    int njetmatched = 0, njetunmatched = 0;
    int index_minjmetdphi = -1, index_maxjmetdphi = -1;
    int index_maxjetptprojonmet_plus = -1, index_maxjetptprojonmet_minus = -1;

    std::vector<int> jetsetaminuseta2p6to3p2;
    std::vector<int> jetsetapluseta2p6to3p2;
    std::vector<int> jetsetaminuseta0to4p7;
    std::vector<int> jetsetapluseta0to4p7;

    for(auto it : *active_part->at(ePos)) {
      histAddVal(part->p4(it).Energy(), "Energy");
      histAddVal(part->p4(it).Pt(), "Pt");
      histAddVal(part->p4(it).Eta(), "Eta");
      histAddVal(part->p4(it).Phi(), "Phi");
      histAddVal(part->p4(it).DeltaPhi(_MET->p4()), "MetDphi");
      if(part->type == PType::Tau) {
        histAddVal(_Tau->charge(it), "Charge");
      }
      if(part->type != PType::Jet) {
        histAddVal(calculateLeptonMetMt(part->p4(it)), "MetMt");
      }
      if(part->type == PType::FatJet ) {
        histAddVal(_FatJet->PrunedMass[it], "PrunedMass");
        histAddVal(_FatJet->SoftDropMass[it], "SoftDropMass");
        histAddVal(_FatJet->tau1[it], "tau1");
        histAddVal(_FatJet->tau2[it], "tau2");
        histAddVal(_FatJet->tau2[it]/_FatJet->tau1[it], "tau2Overtau1");
      }

      if(part->type == PType::Jet){

        // Find out if this is the closest jet to the MET
        float deltaPhiMet = absnormPhi(part->p4(it).Phi() - _MET->phi());
        histAddVal(deltaPhiMet, "AbsDPhiMet");
        histAddVal(normPhi(part->p4(it).Phi() - _MET->phi()), "DPhiMet");
        histAddVal2(normPhi(part->p4(it).Phi() - _MET->phi()), _MET->pt(), "MetVsDPhiMet");

        float cosDPhiJetMet = cos(normPhi(part->p4(it).Phi() - _MET->phi()));
        float jetptmetproj_plus = 0.0, jetptmetproj_minus = 0.0;

        if(cosDPhiJetMet < 0.0){
          jetptmetproj_minus = part->p4(it).Pt() * cosDPhiJetMet;
        } else if(cosDPhiJetMet > 0.0){
          jetptmetproj_plus = part->p4(it).Pt() * cosDPhiJetMet;
        }

        if(deltaPhiMet < minDeltaPhiMet){
          minDeltaPhiMet = deltaPhiMet;
          index_minjmetdphi = it;
        }

        if(deltaPhiMet > maxDeltaPhiMet){
          maxDeltaPhiMet = deltaPhiMet;
          index_maxjmetdphi = it;
        }

        if(abs(jetptmetproj_minus) > maxjetptprojonmet_minus){
          maxjetptprojonmet_minus = abs(jetptmetproj_minus);
          index_maxjetptprojonmet_minus = it;
        }

        if(abs(jetptmetproj_plus) > maxjetptprojonmet_plus){
          maxjetptprojonmet_plus = abs(jetptmetproj_plus);
          index_maxjetptprojonmet_plus = it;
        }

        // ---------
        float jetrawpt = _Jet->RecoP4(it).Pt() * ( 1.0 - _Jet->rawFactor[it] );
        histAddVal2(part->p4(it).Eta(), part->p4(it).Pt(), "PtVsEta");
        histAddVal2(part->p4(it).Eta(), normPhi(part->p4(it).Phi() - _MET->phi()), "MetDphiVsEta");
        histAddVal(jetrawpt, "RawPt");
	       histAddVal2(_Jet->eta(it), jetrawpt, "RawPtvsEta");
        histAddVal2(_Jet->eta(it), _Jet->phi(it), "PhivsEta");
        histAddVal2(_Jet->phi(it), jetrawpt, "RawPtvsPhi");
        histAddVal2(jetrawpt, _Jet->pt(it), "PtvsRawPt");

        /* -----  Histograms to compare neutral/charged EM and hadronic energy fractions */

        float totalHadronicEF = _Jet->neutralHadEnergyFraction[it] + _Jet->chargedHadronEnergyFraction[it];
        float totalEmEF = _Jet->neutralEmEnergyFraction[it] + _Jet->chargedEmEnergyFraction[it];
        float totalneutralEF = _Jet->neutralEmEnergyFraction[it] + _Jet->neutralHadEnergyFraction[it];
        float totalchargedEF = _Jet->chargedEmEnergyFraction[it] + _Jet->chargedHadronEnergyFraction[it];

        if(_Jet->eta(it) > 0.0){
          jetsetapluseta0to4p7.push_back(it);
        } else if (_Jet->eta(it) < 0.0){
          jetsetaminuseta0to4p7.push_back(it);
        }

        if((_Jet->eta(it) > -3.15) && (_Jet->eta(it) < -2.66)){

          jetsetaminuseta2p6to3p2.push_back(it);
          // Reference cold cell
          if( (_Jet->phi(it) > -2.8) && (_Jet->phi(it) < -2.17) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "ColdCellEtaMneHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "ColdCellEtaMneEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "ColdCellEtaMchHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "ColdCellEtaMchEmEF");
            histAddVal(totalHadronicEF                     , "ColdCellEtaMTotalHadEF");
            histAddVal(totalEmEF                           , "ColdCellEtaMTotalEmEF");
            histAddVal(totalneutralEF                      , "ColdCellEtaMTotalNeutralEF");
            histAddVal(totalchargedEF                      , "ColdCellEtaMTotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "ColdCellEtaMnumberOfConstituents");
            histAddVal(jetrawpt                            , "ColdCellEtaMRawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "ColdCellEtaMneHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "ColdCellEtaMneEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "ColdCellEtaMchHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "ColdCellEtaMchEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                      ), "ColdCellEtaMTotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                            ), "ColdCellEtaMTotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                       ), "ColdCellEtaMTotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                       ), "ColdCellEtaMTotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]       ), "ColdCellEtaMnumberOfConstituentsxRawPt");
            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "ColdCellEtaMneHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "ColdCellEtaMneEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "ColdCellEtaMchHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "ColdCellEtaMchEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "ColdCellEtaMTotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "ColdCellEtaMTotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "ColdCellEtaMTotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "ColdCellEtaMTotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "ColdCellEtaMnumberOfConstituentsvsRawPt");
            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt()         , "ColdCellEtaMPtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt                       , "ColdCellEtaMRawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)                  , "ColdCellEtaMPhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt                       , "ColdCellEtaMRawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)                        , "ColdCellEtaMPtvsRawPt");

          }
          // Hot cell, central phi
          if( (_Jet->phi(it) > -0.42) && (_Jet->phi(it) < 0.21) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "HotCellEtaMPhi0neHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "HotCellEtaMPhi0neEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "HotCellEtaMPhi0chHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "HotCellEtaMPhi0chEmEF");
            histAddVal(totalHadronicEF                     , "HotCellEtaMPhi0TotalHadEF");
            histAddVal(totalEmEF                           , "HotCellEtaMPhi0TotalEmEF");
            histAddVal(totalneutralEF                      , "HotCellEtaMPhi0TotalNeutralEF");
            histAddVal(totalchargedEF                      , "HotCellEtaMPhi0TotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "HotCellEtaMPhi0numberOfConstituents");
            histAddVal(jetrawpt                            , "HotCellEtaMPhi0RawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "HotCellEtaMPhi0neHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "HotCellEtaMPhi0neEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "HotCellEtaMPhi0chHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "HotCellEtaMPhi0chEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "HotCellEtaMPhi0TotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "HotCellEtaMPhi0TotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "HotCellEtaMPhi0TotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "HotCellEtaMPhi0TotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "HotCellEtaMPhi0numberOfConstituentsxRawPt");

            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt(), "HotCellEtaMPhi0PtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt              , "HotCellEtaMPhi0RawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)         , "HotCellEtaMPhi0PhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt              , "HotCellEtaMPhi0RawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)               , "HotCellEtaMPhi0PtvsRawPt");

            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "HotCellEtaMPhi0neHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "HotCellEtaMPhi0neEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "HotCellEtaMPhi0chHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "HotCellEtaMPhi0chEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "HotCellEtaMPhi0TotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "HotCellEtaMPhi0TotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "HotCellEtaMPhi0TotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "HotCellEtaMPhi0TotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "HotCellEtaMPhi0numberOfConstituentsvsRawPt");

          }

          // Hot cell, phi plus
          if( (_Jet->phi(it) > 1.19) && (_Jet->phi(it) < 2.03) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "HotCellEtaMPhiPneHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "HotCellEtaMPhiPneEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "HotCellEtaMPhiPchHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "HotCellEtaMPhiPchEmEF");
            histAddVal(totalHadronicEF                     , "HotCellEtaMPhiPTotalHadEF");
            histAddVal(totalEmEF                           , "HotCellEtaMPhiPTotalEmEF");
            histAddVal(jetrawpt                            , "HotCellEtaMPhiPRawPt");
            histAddVal(totalneutralEF                      , "HotCellEtaMPhiPTotalNeutralEF");
            histAddVal(totalchargedEF                      , "HotCellEtaMPhiPTotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "HotCellEtaMPhiPnumberOfConstituents");

            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "HotCellEtaMPhiPneHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "HotCellEtaMPhiPneEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "HotCellEtaMPhiPchHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "HotCellEtaMPhiPchEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "HotCellEtaMPhiPTotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "HotCellEtaMPhiPTotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "HotCellEtaMPhiPTotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "HotCellEtaMPhiPTotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "HotCellEtaMPhiPnumberOfConstituentsxRawPt");


            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt(), "HotCellEtaMPhiPPtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt              , "HotCellEtaMPhiPRawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)         , "HotCellEtaMPhiPPhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt              , "HotCellEtaMPhiPRawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)               , "HotCellEtaMPhiPPtvsRawPt");

            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "HotCellEtaMPhiPneHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "HotCellEtaMPhiPneEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "HotCellEtaMPhiPchHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "HotCellEtaMPhiPchEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                      , "HotCellEtaMPhiPTotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                            , "HotCellEtaMPhiPTotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                       , "HotCellEtaMPhiPTotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                       , "HotCellEtaMPhiPTotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]       , "HotCellEtaMPhiPnumberOfConstituentsvsRawPt");
          }

        }

        if((_Jet->eta(it) > 2.66) && (_Jet->eta(it) < 3.15)){

          jetsetapluseta2p6to3p2.push_back(it);
          // Reference cold cell
          if( (_Jet->phi(it) > -2.8) && (_Jet->phi(it) < -2.17) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "ColdCellEtaPneHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "ColdCellEtaPneEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "ColdCellEtaPchHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "ColdCellEtaPchEmEF");
            histAddVal(totalHadronicEF                     , "ColdCellEtaPTotalHadEF");
            histAddVal(totalEmEF                           , "ColdCellEtaPTotalEmEF");
            histAddVal(totalneutralEF                      , "ColdCellEtaPTotalNeutralEF");
            histAddVal(totalchargedEF                      , "ColdCellEtaPTotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "ColdCellEtaPnumberOfConstituents");
            histAddVal(jetrawpt                            , "ColdCellEtaPRawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "ColdCellEtaPneHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "ColdCellEtaPneEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "ColdCellEtaPchHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "ColdCellEtaPchEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "ColdCellEtaPTotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "ColdCellEtaPTotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "ColdCellEtaPTotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "ColdCellEtaPTotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "ColdCellEtaPnumberOfConstituentsxRawPt");

            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "ColdCellEtaPneHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "ColdCellEtaPneEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "ColdCellEtaPchHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "ColdCellEtaPchEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "ColdCellEtaPTotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "ColdCellEtaPTotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "ColdCellEtaPTotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "ColdCellEtaPTotalChargedEFvsRawPt");
            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt()         , "ColdCellEtaPPtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt                       , "ColdCellEtaPRawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)                  , "ColdCellEtaPPhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt                       , "ColdCellEtaPRawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)                        , "ColdCellEtaPPtvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "ColdCellEtaPnumberOfConstituentsvsRawPt");

          }

          // Hot cell, phi minus
          if( (_Jet->phi(it) > -2.03) && (_Jet->phi(it) < -1.19) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhiMneHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhiMneEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhiMchHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhiMchEmEF");
            histAddVal(totalHadronicEF                     , "HotCellEtaPPhiMTotalHadEF");
            histAddVal(totalEmEF                           , "HotCellEtaPPhiMTotalEmEF");
            histAddVal(totalneutralEF                      , "HotCellEtaPPhiMTotalNeutralEF");
            histAddVal(totalchargedEF                      , "HotCellEtaPPhiMTotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "HotCellEtaPPhiMnumberOfConstituents");

            histAddVal(jetrawpt                            , "HotCellEtaPPhiMRawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "HotCellEtaPPhiMneHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "HotCellEtaPPhiMneEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "HotCellEtaPPhiMchHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "HotCellEtaPPhiMchEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "HotCellEtaPPhiMTotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "HotCellEtaPPhiMTotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "HotCellEtaPPhiMTotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "HotCellEtaPPhiMTotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "HotCellEtaPPhiMnumberOfConstituentsxRawPt");


            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt(), "HotCellEtaPPhiMPtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt              , "HotCellEtaPPhiMRawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)         , "HotCellEtaPPhiMPhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt              , "HotCellEtaPPhiMRawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)               , "HotCellEtaPPhiMPtvsRawPt");

            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhiMneHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhiMneEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhiMchHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhiMchEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "HotCellEtaPPhiMTotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "HotCellEtaPPhiMTotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "HotCellEtaPPhiMTotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "HotCellEtaPPhiMTotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "HotCellEtaPPhiMnumberOfConstituentsvsRawPt");

          }


          // Hot cell, central phi
          if( (_Jet->phi(it) > -0.42) && (_Jet->phi(it) < 0.42) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhi0neHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhi0neEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhi0chHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhi0chEmEF");
            histAddVal(totalHadronicEF                     , "HotCellEtaPPhi0TotalHadEF");
            histAddVal(totalEmEF                           , "HotCellEtaPPhi0TotalEmEF");
            histAddVal(totalneutralEF                      , "HotCellEtaPPhi0TotalNeutralEF");
            histAddVal(totalchargedEF                      , "HotCellEtaPPhi0TotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "HotCellEtaPPhi0numberOfConstituents");
            histAddVal(jetrawpt                            , "HotCellEtaPPhi0RawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "HotCellEtaPPhi0neHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "HotCellEtaPPhi0neEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "HotCellEtaPPhi0chHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "HotCellEtaPPhi0chEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "HotCellEtaPPhi0TotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "HotCellEtaPPhi0TotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "HotCellEtaPPhi0TotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "HotCellEtaPPhi0TotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "HotCellEtaPPhi0numberOfConstituentsxRawPt");

            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt(), "HotCellEtaPPhi0PtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt              , "HotCellEtaPPhi0RawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)         , "HotCellEtaPPhi0PhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt              , "HotCellEtaPPhi0RawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)               , "HotCellEtaPPhi0PtvsRawPt");

            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhi0neHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhi0neEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhi0chHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhi0chEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "HotCellEtaPPhi0TotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "HotCellEtaPPhi0TotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "HotCellEtaPPhi0TotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "HotCellEtaPPhi0TotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "HotCellEtaPPhi0numberOfConstituentsvsRawPt");
          }

          // Hot cell, phi plus
          if( (_Jet->phi(it) > 1.05) && (_Jet->phi(it) < 1.61) ){

            histAddVal(_Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhiPneHEF");
            histAddVal(_Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhiPneEmEF");
            histAddVal(_Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhiPchHEF");
            histAddVal(_Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhiPchEmEF");
            histAddVal(totalHadronicEF                     , "HotCellEtaPPhiPTotalHadEF");
            histAddVal(totalEmEF                           , "HotCellEtaPPhiPTotalEmEF");
            histAddVal(totalneutralEF                      , "HotCellEtaPPhiPTotalNeutralEF");
            histAddVal(totalchargedEF                      , "HotCellEtaPPhiPTotalChargedEF");
            histAddVal(_Jet->numberOfConstituents[it]      , "HotCellEtaPPhiPnumberOfConstituents");
            histAddVal(jetrawpt                            , "HotCellEtaPPhiPRawPt");
            histAddVal( (jetrawpt * _Jet->neutralHadEnergyFraction[it]   ), "HotCellEtaPPhiPneHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->neutralEmEnergyFraction[it]    ), "HotCellEtaPPhiPneEmEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedHadronEnergyFraction[it]), "HotCellEtaPPhiPchHEFxRawPt");
            histAddVal( (jetrawpt * _Jet->chargedEmEnergyFraction[it]    ), "HotCellEtaPPhiPchEmEFxRawPt");
            histAddVal( (jetrawpt * totalHadronicEF                     ), "HotCellEtaPPhiPTotalHadEFxRawPt");
            histAddVal( (jetrawpt * totalEmEF                           ), "HotCellEtaPPhiPTotalEmEFxRawPt");
            histAddVal( (jetrawpt * totalneutralEF                      ), "HotCellEtaPPhiPTotalNeutralEFxRawPt");
            histAddVal( (jetrawpt * totalchargedEF                      ), "HotCellEtaPPhiPTotalChargedEFxRawPt");
            histAddVal( (jetrawpt * _Jet->numberOfConstituents[it]      ), "HotCellEtaPPhiPnumberOfConstituentsxRawPt");

            histAddVal2(_Jet->p4(it).Eta(), _Jet->p4(it).Pt(), "HotCellEtaPPhiPPtVsEta");
            histAddVal2(_Jet->eta(it), jetrawpt              , "HotCellEtaPPhiPRawPtvsEta");
            histAddVal2(_Jet->eta(it), _Jet->phi(it)         , "HotCellEtaPPhiPPhivsEta");
            histAddVal2(_Jet->phi(it), jetrawpt              , "HotCellEtaPPhiPRawPtvsPhi");
            histAddVal2(jetrawpt, _Jet->pt(it)               , "HotCellEtaPPhiPPtvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralHadEnergyFraction[it]   , "HotCellEtaPPhiPneHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->neutralEmEnergyFraction[it]    , "HotCellEtaPPhiPneEmEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedHadronEnergyFraction[it], "HotCellEtaPPhiPchHEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->chargedEmEnergyFraction[it]    , "HotCellEtaPPhiPchEmEFvsRawPt");
            histAddVal2(jetrawpt, totalHadronicEF                     , "HotCellEtaPPhiPTotalHadEFvsRawPt");
            histAddVal2(jetrawpt, totalEmEF                           , "HotCellEtaPPhiPTotalEmEFvsRawPt");
            histAddVal2(jetrawpt, totalneutralEF                      , "HotCellEtaPPhiPTotalNeutralEFvsRawPt");
            histAddVal2(jetrawpt, totalchargedEF                      , "HotCellEtaPPhiPTotalChargedEFvsRawPt");
            histAddVal2(jetrawpt, _Jet->numberOfConstituents[it]      , "HotCellEtaPPhiPnumberOfConstituentsvsRawPt");
          }

        }

        /* ----------------------------------------------------------------------------  */


        if( _Jet->pstats["Smear"].bfind("SmearTheJet") ){
          std::map<int, std::vector<float> >::iterator smdj = jets_jer_sfs.find(it);

          if(smdj != jets_jer_sfs.end()){
            if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 0 && _Jet->getPileupJetID(it,2) == 1){
              // Pass Loose PU jet ID
              histAddVal(3, "PUJetID");
            } else if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 1 ){
              // Pass Medium PU jet ID
              histAddVal(2, "PUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 1){
              // Pass Tight PU jet ID
              histAddVal(1, "PUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0){
              // Fail PU jet ID
              histAddVal(0, "PUJetID");
            }

            // Histograms for general purpose:
            // Nominal JER SF -- c_JER: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
            histAddVal2(part->p4(it).Eta(), jets_jer_sfs[it].at(0), "cJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jer_sfs[it].at(0), "cJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jer_sfs[it].at(0), "cJERNomvsPhi");
            // Nominal data to simulation core resolution scale factor, s_JER
            histAddVal2(part->p4(it).Eta(), jets_datatoMCSFs_jer[it], "sJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_datatoMCSFs_jer[it], "sJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_datatoMCSFs_jer[it], "sJERNomvsPhi");
            // Nominal relative pt resolution, sigma_JER
            histAddVal2(part->p4(it).Eta(), jets_jetptres[it], "sigmaJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jetptres[it], "sigmaJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jetptres[it], "sigmaJERNomvsPhi");
          }

          std::vector<int>::iterator gmj = std::find(genMatchedJets.begin(), genMatchedJets.end(), it);
          std::vector<int>::iterator gumj = std::find(genUnmatchedJets.begin(), genUnmatchedJets.end(), it);


          if(gmj != genMatchedJets.end()){

            histAddVal(jets_jer_sfs[it].at(0), "GenMatchedcJERNom");
            histAddVal(jets_datatoMCSFs_jer[it], "GenMatchedsJERNom");
            histAddVal(jets_jetptres[it], "GenMatchedsigmaJERNom");
            histAddVal(part->p4(it).Pt(), "GenMatchedPt");
            histAddVal(part->p4(it).Eta(), "GenMatchedEta");
            histAddVal(part->p4(it).Phi(), "GenMatchedPhi");
            histAddVal(normPhi(part->p4(it).Phi() - _MET->phi()), "GenMatchedMetDphi");

            if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 0 && _Jet->getPileupJetID(it,2) == 1){
              // Pass Loose PU jet ID
              histAddVal(3, "GenMatchedPUJetID");
            } else if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 1 ){
              // Pass Medium PU jet ID
              histAddVal(2, "GenMatchedPUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 1){
              // Pass Tight PU jet ID
              histAddVal(1, "GenMatchedPUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0){
              // Fail PU jet ID
              histAddVal(0, "GenMatchedPUJetID");
            }

            histAddVal2(part->p4(it).Eta(), part->p4(it).Pt(), "GenMatchedPtVsEta");
            histAddVal2(part->p4(it).Eta(), normPhi(part->p4(it).Phi() - _MET->phi()), "GenMatchedMetDphiVsEta");

            // Nominal JER SF -- c_JER:
            histAddVal2(part->p4(it).Eta(), jets_jer_sfs[it].at(0), "GenMatchedcJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jer_sfs[it].at(0), "GenMatchedcJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jer_sfs[it].at(0), "GenMatchedcJERNomvsPhi");
            // Nominal data to simulation core resolution scale factor, s_JER
            histAddVal2(part->p4(it).Eta(), jets_datatoMCSFs_jer[it], "GenMatchedsJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_datatoMCSFs_jer[it], "GenMatchedsJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_datatoMCSFs_jer[it], "GenMatchedsJERNomvsPhi");
            // Nominal relative pt resolution, sigma_JER
            histAddVal2(part->p4(it).Eta(), jets_jetptres[it], "GenMatchedsigmaJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jetptres[it], "GenMatchedsigmaJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jetptres[it], "GenMatchedsigmaJERNomvsPhi");

            histAddVal(deltaPhiMet, "GenMatchedAbsDPhiMet");
            histAddVal(normPhi(part->p4(it).Phi() - _MET->phi()), "GenMatchedDPhiMet");
            histAddVal2(normPhi(part->p4(it).Phi() - _MET->phi()), _MET->pt(), "GenMatchedMetVsDPhiMet");

            njetmatched++;
          } else if(gumj != genUnmatchedJets.end()){
            histAddVal(jets_jer_sfs[it].at(0), "GenUnmatchedcJERNom");
            histAddVal(jets_datatoMCSFs_jer[it], "GenUnmatchedsJERNom");
            histAddVal(jets_jetptres[it], "GenUnmatchedsigmaJERNom");
            histAddVal(part->p4(it).Pt(), "GenUnmatchedPt");
            histAddVal(part->p4(it).Eta(), "GenUnmatchedEta");
            histAddVal(part->p4(it).Phi(), "GenUnmatchedPhi");
            histAddVal(normPhi(part->p4(it).Phi() - _MET->phi()), "GenUnmatchedMetDphi");

            if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 0 && _Jet->getPileupJetID(it,2) == 1){
              // Pass Loose PU jet ID
              histAddVal(3, "GenUnmatchedPUJetID");
            } else if( _Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0 && _Jet->getPileupJetID(it,1) == 1 ){
              // Pass Medium PU jet ID
              histAddVal(2, "GenUnmatchedPUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 1){
              // Pass Tight PU jet ID
              histAddVal(1, "GenUnmatchedPUJetID");
            } else if(_Jet->pt(it) <= 50.0 && _Jet->getPileupJetID(it,0) == 0){
              // Fail PU jet ID
              histAddVal(0, "GenUnmatchedPUJetID");
            }

            histAddVal2(part->p4(it).Eta(), part->p4(it).Pt(), "GenUnmatchedPtVsEta");
            histAddVal2(part->p4(it).Eta(), normPhi(part->p4(it).Phi() - _MET->phi()), "GenUnmatchedMetDphiVsEta");
            // Nominal JER SF -- c_JER:
            histAddVal2(part->p4(it).Eta(), jets_jer_sfs[it].at(0), "GenUnmatchedcJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jer_sfs[it].at(0), "GenUnmatchedcJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jer_sfs[it].at(0), "GenUnmatchedcJERNomvsPhi");
            // Nominal data to simulation core resolution scale factor, s_JER
            histAddVal2(part->p4(it).Eta(), jets_datatoMCSFs_jer[it], "GenUnmatchedsJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_datatoMCSFs_jer[it], "GenUnmatchedsJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_datatoMCSFs_jer[it], "GenUnmatchedsJERNomvsPhi");
            // Nominal relative pt resolution, sigma_JER
            histAddVal2(part->p4(it).Eta(), jets_jetptres[it], "GenUnmatchedsigmaJERNomvsEta");
            histAddVal2(part->p4(it).Pt(), jets_jetptres[it], "GenUnmatchedsigmaJERNomvsPt");
            histAddVal2(part->p4(it).Phi(), jets_jetptres[it], "GenUnmatchedsigmaJERNomvsPhi");

            histAddVal(deltaPhiMet, "GenUnmatchedAbsDPhiMet");
            histAddVal(normPhi(part->p4(it).Phi() - _MET->phi()), "GenUnmatchedDPhiMet");
            histAddVal2(normPhi(part->p4(it).Phi() - _MET->phi()), _MET->pt(), "GenUnmatchedMetVsDPhiMet");

            njetunmatched++;
          }


        }
        // ---------
      }

      if(ePos == CUTS::eRBJet){
        histAddVal2(part->p4(it).Eta(), part->p4(it).Pt(), "PtVsEta");
      }

    }

    if(part->type == PType::Jet){
      histAddVal(minDeltaPhiMet, "MinAbsDPhiMet"); // minimum |deltaPhi(jet, MET)|

      if(index_minjmetdphi > -1){
        histAddVal(_Jet->pt(index_minjmetdphi), "MinAbsDPhiMetPt");
        histAddVal(_Jet->eta(index_minjmetdphi), "MinAbsDPhiMetEta");
        histAddVal2(_Jet->eta(index_minjmetdphi), _MET->pt(), "MetvsMinAbsDPhiMetEta");
        histAddVal2(_Jet->eta(index_minjmetdphi), _Jet->pt(index_minjmetdphi), "MetvsMinAbsDPhiMetPtvsEta");
        histAddVal(_Jet->phi(index_minjmetdphi), "MinAbsDPhiMetPhi");
        histAddVal(minDeltaPhiMet, "MinAbsDPhiMet");
        histAddVal2(minDeltaPhiMet, _MET->pt(), "MetvsMinAbsDPhiMet");
      }

      if(index_maxjmetdphi > -1){
        histAddVal(_Jet->pt(index_maxjmetdphi), "MaxAbsDPhiMetPt");
        histAddVal(_Jet->eta(index_maxjmetdphi), "MaxAbsDPhiMetEta");
        histAddVal2(_Jet->eta(index_maxjmetdphi), _MET->pt(), "MetvsMaxAbsDPhiMetEta");
        histAddVal2(_Jet->eta(index_maxjmetdphi), _Jet->pt(index_maxjmetdphi), "MetvsMaxAbsDPhiMetPtvsEta");
        histAddVal(_Jet->phi(index_maxjmetdphi), "MaxAbsDPhiMetPhi");
        histAddVal(maxDeltaPhiMet, "MaxAbsDPhiMet");
        histAddVal2(maxDeltaPhiMet, _MET->pt(), "MetvsMaxAbsDPhiMet");
      }

      if(index_maxjetptprojonmet_minus > -1){
        histAddVal(_Jet->pt(index_maxjetptprojonmet_minus), "MaxProjMetMinusPt");
        histAddVal(_Jet->eta(index_maxjetptprojonmet_minus), "MaxProjMetMinusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_minus), _MET->pt(), "MetvsMaxProjMetMinusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_minus), _Jet->pt(index_maxjetptprojonmet_minus), "MetvsMaxProjMetMinusPtvsEta");
        histAddVal(_Jet->phi(index_maxjetptprojonmet_minus), "MaxProjMetMinusPhi");
        histAddVal(maxjetptprojonmet_minus, "MaxProjMetMinusMag");
        histAddVal2(maxjetptprojonmet_minus, _MET->pt(), "MetvsMaxProjMetMinusMag");
      }

      if(index_maxjetptprojonmet_plus > -1){
        histAddVal(_Jet->pt(index_maxjetptprojonmet_plus), "MaxProjMetPlusPt");
        histAddVal(_Jet->eta(index_maxjetptprojonmet_plus), "MaxProjMetPlusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_plus), _MET->pt(), "MetvsMaxProjMetPlusEta");
        histAddVal2(_Jet->eta(index_maxjetptprojonmet_plus), _Jet->pt(index_maxjetptprojonmet_plus), "MetvsMaxProjMetPlusPtvsEta");
        histAddVal(_Jet->phi(index_maxjetptprojonmet_plus), "MaxProjMetPlusPhi");
        histAddVal(maxjetptprojonmet_plus, "MaxProjMetPlusMag");
        histAddVal2(maxjetptprojonmet_plus, _MET->pt(), "MetvsMaxProjMetPlusMag");
      }

      if( index_maxjetptprojonmet_minus > -1 || index_maxjetptprojonmet_plus > -1){
        int maxprojonmet_idx = maxjetptprojonmet_minus > maxjetptprojonmet_plus ? index_maxjetptprojonmet_minus : index_maxjetptprojonmet_plus;
        float maxprojonmet = std::max(maxjetptprojonmet_minus, maxjetptprojonmet_plus);

        histAddVal(_Jet->pt(maxprojonmet_idx), "MaxProjMetPt");
        histAddVal(_Jet->eta(maxprojonmet_idx), "MaxProjMetEta");
        histAddVal2(_Jet->eta(maxprojonmet_idx), _MET->pt(), "MetvsMaxProjMetEta");
        histAddVal2(_Jet->eta(maxprojonmet_idx), _Jet->pt(maxprojonmet_idx), "MetvsMaxProjMetPtvsEta");
        histAddVal(_Jet->phi(maxprojonmet_idx), "MaxProjMetPhi");
        histAddVal(maxprojonmet, "MaxProjMetMag");
        histAddVal2(maxprojonmet, _MET->pt(), "MetvsMaxProjMetMag");
        histAddVal(normPhi(_Jet->phi(maxprojonmet_idx) - _MET->phi()), "MaxProjMetDeltaPhi");
        histAddVal2(normPhi(_Jet->phi(maxprojonmet_idx) - _MET->phi()), _Jet->pt(maxprojonmet_idx), "MetvsMaxProjMetDeltaPhi");
        if( cos( normPhi(_Jet->phi(maxprojonmet_idx) - _MET->phi()) ) > 0 ){
          histAddVal(1, "MaxProjMetSign");
          histAddVal2(1, _MET->pt(), "MetvsMaxProjMetSign");
        } else {
          histAddVal(-1, "MaxProjMetSign");
          histAddVal2(-1, _MET->pt(), "MetvsMaxProjMetSign");
        }
        histAddVal(cos( normPhi(_Jet->phi(maxprojonmet_idx) - _MET->phi()) ), "MaxProjMetCosDPhi");
      }

      histAddVal2(minDeltaPhiMet, _MET->pt(), "MetVsMinAbsDPhiMet");
      histAddVal(njetmatched, "NGenMatched");
      histAddVal(njetunmatched, "NGenUnmatched");
      histAddVal2(active_part->at(ePos)->size(), njetmatched, "NGenMatchedVsN");
      histAddVal2(active_part->at(ePos)->size(), njetunmatched, "NGenMatchedVsN");
      histAddVal2(njetunmatched, njetmatched, "NGenMatchedVsNGenUnmatched");
      histAddVal2(njetmatched, njetunmatched, "NGenUnmatchedVsNGenMatched");
      histAddVal2(_MET->pt(), njetmatched, "NGenMatchedVsMetPt");
      histAddVal2(_MET->pt(), njetunmatched, "NGenUnmatchedVsMetPt");
      histAddVal2(_MET->phi(), njetmatched, "NGenMatchedVsMetPhi");
      histAddVal2(_MET->phi(), njetunmatched, "NGenUnmatchedVsMetPhi");
      histAddVal2(bestVertices, njetmatched, "NGenMatchedVsNPV");
      histAddVal2(bestVertices, njetunmatched, "NGenUnmatchedVsNPV");

      std::vector<int>::iterator gmj = std::find(genMatchedJets.begin(), genMatchedJets.end(), index_minjmetdphi);
      std::vector<int>::iterator gumj = std::find(genUnmatchedJets.begin(), genUnmatchedJets.end(), index_minjmetdphi);

      if(gmj != genMatchedJets.end()){
        histAddVal(minDeltaPhiMet, "GenMatchedMinAbsDPhiMet");
        histAddVal2(minDeltaPhiMet, _MET->pt(), "GenMatchedMetVsMinAbsDPhiMet");
      } else if(gumj != genUnmatchedJets.end()){
        histAddVal(minDeltaPhiMet, "GenUnmatchedMinAbsDPhiMet");
        histAddVal2(minDeltaPhiMet, _MET->pt(), "GenUnmatchedMetVsMinAbsDPhiMet");

      }

      histAddVal(jetsetapluseta2p6to3p2.size(), "NPosEEnoiseEta");
      histAddVal(jetsetaminuseta2p6to3p2.size(), "NNegEEnoiseEta");
      
      histAddVal2(jetsetaminuseta0to4p7.size(), jetsetapluseta2p6to3p2.size(), "NEEnoiseEtaPosvsTotalNeg");
      histAddVal2(jetsetaminuseta2p6to3p2.size(), jetsetapluseta0to4p7.size(), "NEEnoiseEtaTotalPosvsNeg");
      histAddVal2(jetsetaminuseta2p6to3p2.size(), jetsetapluseta2p6to3p2.size(), "NEEnoiseEtaPosvsNeg");
      histAddVal2(jetsetaminuseta0to4p7.size(), jetsetapluseta0to4p7.size(), "NEtaTotalPosvsTotalNeg");

      if(jetsetaminuseta2p6to3p2.size() > 0){
        histAddVal(jetsetapluseta0to4p7.size(), "NNegEEnoiseEtaPlusEta");
      }
      if(jetsetapluseta2p6to3p2.size() > 0){
        histAddVal(jetsetaminuseta0to4p7.size(), "NPlusEEnoiseEtaNegEta");
      }

      if(jetsetapluseta2p6to3p2.size() > 0 && jetsetaminuseta2p6to3p2.size() > 0){
        histAddVal(1, "NOSEEnoiseEta");

        for(size_t ip = 0; ip < jetsetapluseta2p6to3p2.size(); ip++){
          int index_posj = jetsetapluseta2p6to3p2.at(ip);
          
          for(size_t im = 0; im < jetsetaminuseta2p6to3p2.size(); im++){
            int index_negj = jetsetaminuseta2p6to3p2.at(im);

            TLorentzVector deltaPEENoiseJets = (_Jet->p4(index_posj) - _Jet->p4(index_negj));

            float totalEmEF  = _Jet->neutralEmEnergyFraction[index_posj] + _Jet->chargedEmEnergyFraction[index_posj] + _Jet->neutralEmEnergyFraction[index_negj] + _Jet->chargedEmEnergyFraction[index_negj];
            float totalHadEF = _Jet->neutralHadEnergyFraction[index_posj] + _Jet->chargedHadronEnergyFraction[index_posj] + _Jet->neutralHadEnergyFraction[index_negj] + _Jet->chargedHadronEnergyFraction[index_negj];

            float totalneutralEmEF  = (_Jet->neutralEmEnergyFraction[index_posj] + _Jet->neutralEmEnergyFraction[index_negj]) / totalEmEF;
            float totalchargedEmEF  = (_Jet->chargedEmEnergyFraction[index_posj] + _Jet->chargedEmEnergyFraction[index_negj]) / totalEmEF;
            float totalneutralHadEF = (_Jet->neutralHadEnergyFraction[index_posj] + _Jet->neutralHadEnergyFraction[index_negj]) / totalHadEF;
            float totalchargedHadEF = (_Jet->chargedHadronEnergyFraction[index_posj] + _Jet->chargedHadronEnergyFraction[index_negj]) / totalHadEF;

            histAddVal(deltaPEENoiseJets.Pt()                    , "EEnoiseCombDeltaPt");
            histAddVal((deltaPEENoiseJets.Pt() / _MET->HT())     , "EEnoiseCombDeltaPtHTRatio");
            histAddVal2(_MET->HT(), deltaPEENoiseJets.Pt()       , "EEnoiseCombDeltaPtvsHT");
            histAddVal((deltaPEENoiseJets.Pt() / _MET->pt())     , "EEnoiseCombDeltaPtMetRatio");
            histAddVal2(_MET->pt(), deltaPEENoiseJets.Pt()       , "EEnoiseCombDeltaPtvsMet");
            histAddVal2(deltaPEENoiseJets.Pt(), totalEmEF        , "EEnoiseCombTotalEmEFvsDeltaPt");
            histAddVal2(deltaPEENoiseJets.Pt(), totalHadEF       , "EEnoiseCombTotalHadEFvsDeltaPt");
            histAddVal2(deltaPEENoiseJets.Pt(), totalneutralEmEF , "EEnoiseCombTotalNeutralEmEFvsDeltaPt");
            histAddVal2(deltaPEENoiseJets.Pt(), totalchargedEmEF , "EEnoiseCombTotalChargedEmEFvsDeltaPt");
            histAddVal2(deltaPEENoiseJets.Pt(), totalneutralHadEF, "EEnoiseCombTotalNeutralHadEFvsDeltaPt");
            histAddVal2(deltaPEENoiseJets.Pt(), totalchargedHadEF, "EEnoiseCombTotalChargedHadEFvsDeltaPt");

          } 
        }
      } else{
        histAddVal(0, "NOSEEnoiseEta");
      }
      
    }

    if((part->type != PType::Jet ) && active_part->at(ePos)->size() > 0) {
      std::vector<std::pair<double, int> > ptIndexVector;
      for(auto it : *active_part->at(ePos)) {
        ptIndexVector.push_back(std::make_pair(part->pt(it),it));
      }
      sort(ptIndexVector.begin(),ptIndexVector.end());
      if(ptIndexVector.size()>0){
        histAddVal(part->pt(ptIndexVector.back().second), "FirstLeadingPt");
        histAddVal(part->eta(ptIndexVector.back().second), "FirstLeadingEta");
      	if((part->type == PType::Muon )){
      		Float_t  leadingmu_pt = (part->pt(ptIndexVector.back().second));
      		Float_t hnandfatjets = leadingmu_pt;
      		for(size_t i=0; i<active_part->at(CUTS::eRWjet)->size(); i++){
      			hnandfatjets = hnandfatjets + (_FatJet->p4(active_part->at(CUTS::eRWjet)->at(i)).Pt());
      		}
      		histAddVal(hnandfatjets, "ptak8pt");
      	}
      }
      if(ptIndexVector.size()>1){
        histAddVal(part->pt(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingPt");
        histAddVal(part->eta(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingEta");
      }
    }

    histAddVal(active_part->at(ePos)->size(), "N");

  } else if(group == "FillMetCuts") {
    histAddVal(_MET->MHT(), "MHT");
    histAddVal(_MET->HT(), "HT");
    histAddVal(_MET->HT() + _MET->MHT(), "Meff");
    histAddVal(_MET->pt(), "Met");
    histAddVal(_MET->phi(), "MetPhi");
    // Added by Brenda - 02/22/2021
    histAddVal2(_MET->pt(), totalVertices, "NPVvsMetPt");
    histAddVal2(_MET->phi(), totalVertices, "NPVvsMetPhi");
    histAddVal(_MET->DefMet.Pt(), "DefaultMETOriginal");
    histAddVal(_MET->T1Met.Pt(), "T1METOriginal");
    histAddVal(_MET->RawMet.Pt(), "RawMETOriginal");
    histAddVal(_MET->JERCorrMet.Pt(),"CorrectedRawMET");

  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() == 0) {

    if(active_part->at(CUTS::eR1stJet)->size()>0) { //01.17.19
      // std::cout << "First leading jet index (fill) = " << active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1) << ", pt = " << _Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Pt() << std::endl;
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Pt(), "FirstPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Eta(), "FirstEta");
      Double_t dphi1new = normPhi(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Phi() - _MET->phi());
      histAddVal(dphi1new,"Dphi1");
    }
    if(active_part->at(CUTS::eR2ndJet)->size()>0) {
      // std::cout << "Second leading jet index (fill) = " << active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1) << ", pt = " << _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1)).Pt() << std::endl;
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1)).Pt(), "SecondPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1)).Eta(), "SecondEta");
    }

  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() != 0) {

    TLorentzVector first = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size() - 1));
    TLorentzVector second = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size() - 1));
    float firstjetrawpt = first.Pt() * (1.0 -  _Jet->rawFactor[active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size() - 1)]);
    float secondjetrawpt = second.Pt() * (1.0 - _Jet->rawFactor[active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size() - 1)] );

    histAddVal(first.Pt(), "FirstPt");
    histAddVal(second.Pt(), "SecondPt");

    histAddVal(first.Eta(), "FirstEta");
    histAddVal(second.Eta(), "SecondEta");

    histAddVal(first.Phi(), "FirstPhi");
    histAddVal(second.Phi(), "SecondPhi");

    histAddVal(firstjetrawpt, "FirstRawPt");
    histAddVal(secondjetrawpt, "SecondRawPt");

    histAddVal2(first.Eta(), first.Pt(), "FirstPtvsEta");
    histAddVal2(first.Phi(), first.Pt(), "FirstPtvsPhi");
    histAddVal2(first.Eta(), first.Phi(), "FirstPhivsEta");
    histAddVal2(firstjetrawpt, first.Pt(), "FirstPtvsRawPt");
    histAddVal2(first.Eta(), firstjetrawpt, "FirstRawPtvsEta");
    histAddVal2(first.Phi(), firstjetrawpt, "FirstRawPtvsPhi");

    histAddVal2(second.Eta(), second.Pt(), "SecondPtvsEta");
    histAddVal2(second.Phi(), second.Pt(), "SecondPtvsPhi");
    histAddVal2(second.Eta(), second.Phi(), "SecondPhivsEta");
    histAddVal2(secondjetrawpt, second.Pt(), "SecondPtvsRawPt");
    histAddVal2(second.Eta(), secondjetrawpt, "SecondRawPtvsEta");
    histAddVal2(second.Phi(), secondjetrawpt, "SecondRawPtvsPhi");


    TLorentzVector LeadDiJet = first + second;

    histAddVal(LeadDiJet.M(), "Mass");
    histAddVal(LeadDiJet.Pt(), "Pt");
    histAddVal(fabs(first.Eta() - second.Eta()), "DeltaEta");
    histAddVal(first.DeltaR(second), "DeltaR");

    float dphiDijets = absnormPhi(first.Phi() - second.Phi());
    float dphi1 = normPhi(first.Phi() - _MET->phi());
    float dphi2 = normPhi(second.Phi() - _MET->phi());
    float alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi");
    histAddVal2(_MET->pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(absnormPhi(_MET->phi() - LeadDiJet.Phi()), "MetDeltaPhi");

    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - _MET->MHTphi()), "Dphi1MHT");
    histAddVal(normPhi(second.Phi() - _MET->MHTphi()), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");

    int idx_flj = active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size() - 1);
    int idx_slj = active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size() - 1);

    if( _Jet->pstats["Smear"].bfind("SmearTheJet") ){
      std::vector<int>::iterator gmj1 = std::find(genMatchedJets.begin(), genMatchedJets.end(), idx_flj);
      std::vector<int>::iterator gumj1 = std::find(genUnmatchedJets.begin(), genUnmatchedJets.end(), idx_flj);

      if(gmj1 != genMatchedJets.end()){
        histAddVal(1, "FirstGenMatching");
        histAddVal(dphi1, "FirstGenMatchedMetDphi");

      } else if(gumj1 != genUnmatchedJets.end()){
        histAddVal(0, "FirstGenMatching");
        histAddVal(dphi1, "FirstGenUnmatchedMetDphi");
      }

      std::vector<int>::iterator gmj2 = std::find(genMatchedJets.begin(), genMatchedJets.end(), idx_slj);
      std::vector<int>::iterator gumj2 = std::find(genUnmatchedJets.begin(), genUnmatchedJets.end(), idx_slj);

      if(gmj2 != genMatchedJets.end()){
        histAddVal(1, "SecondGenMatching");
        histAddVal(dphi2, "SecondGenMatchedMetDphi");

      } else if(gumj2 != genUnmatchedJets.end()){
        histAddVal(0, "SecondGenMatching");
        histAddVal(dphi2, "SecondGenUnmatchedMetDphi");
      }

    }

  } else if(group == "FillDiJet"){ // Dijet combinations
    float leaddijetmass = 0;
    float leaddijetpt = 0;
    float leaddijetdeltaR = 0;
    float leaddijetdeltaEta = 0;
    float etaproduct = 0;
    float largestMassDeltaEta = 0; // added by Kyungmin

    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int p1 = (it) / _Jet->size();
      int p2 = (it) % _Jet->size();
      TLorentzVector jet1 = _Jet->p4(p1);
      TLorentzVector jet2 = _Jet->p4(p2);
      TLorentzVector DiJet = jet1 + jet2;

      if(DiJet.M() > leaddijetmass) {
        leaddijetmass = DiJet.M();
        etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
        largestMassDeltaEta = fabs(jet1.Eta() - jet2.Eta());
      }
      if(DiJet.Pt() > leaddijetpt) leaddijetpt = DiJet.Pt();
      if(fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) leaddijetdeltaEta = fabs(jet1.Eta() - jet2.Eta());
      if(jet1.DeltaR(jet2) > leaddijetdeltaR) leaddijetdeltaR = jet1.DeltaR(jet2);

      histAddVal(DiJet.M(), "Mass");
      histAddVal(DiJet.Pt(), "Pt");
      histAddVal(fabs(jet1.Eta() - jet2.Eta()), "DeltaEta");
      histAddVal(absnormPhi(jet1.Phi() - jet2.Phi()), "DeltaPhi");
      histAddVal(jet1.DeltaR(jet2), "DeltaR");
    }

    histAddVal(leaddijetmass, "LargestMass");
    histAddVal(leaddijetpt, "LargestPt");
    histAddVal(leaddijetdeltaEta, "LargestDeltaEta");
    histAddVal(leaddijetdeltaR, "LargestDeltaR");
    histAddVal(etaproduct, "LargestMassEtaProduct");
    histAddVal(largestMassDeltaEta, "LargestMassDeltaEta");

    for(auto index : *(active_part->at(CUTS::eRTau1)) ) {
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetmass, "mTvsLeadingMass");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaEta, "mTvsLeadingDeltaEta");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaR, "mTvsLeadingDeltaR");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetpt, "mTvsLeadingPt");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetmass, "MetDphiVSLeadingMass");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaEta, "MetDphiVSLeadingDeltaEta");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaR, "MetDphiVSLeadingDeltaR");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetpt, "MetDphiVSLeadingPt");
    }

  } else if(fillInfo[group]->type == FILLER::Dilepjet) { // LeptonJet combinations
    Jet* jet = static_cast<Jet*>(fillInfo[group]->part);
    Lepton* lep = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1(0,0,0,0), part2(0,0,0,0);

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / _Jet->size();;
      int p2= (it) % _Jet->size();;

      part1 = lep->p4(p1);
      part2 = jet->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");
      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }
      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt");
      histAddVal(PZeta, "PZeta");
      histAddVal(PZetaVis, "PZetaVis");
      histAddVal2(PZetaVis,PZeta, "Zeta2D");
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }
    }
  } else if(fillInfo[group]->type == FILLER::Dipart) {   // Dilepton combinations
    Lepton* lep1 = static_cast<Lepton*>(fillInfo[group]->part);
    Lepton* lep2 = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1(0,0,0,0), part2(0,0,0,0), llep(0,0,0,0);

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / BIG_NUM;
      int p2= (it) % BIG_NUM;

      part1 = lep1->p4(p1);
      part2 = lep2->p4(p2);

      // Get their corresponding rest masses
      float mass1 = leptonmasses.at(lep1->type);
      float mass2 = leptonmasses.at(lep2->type);

      if(part1.Pt() > part2.Pt()) llep = lep1->p4(p1);
      else if(part1.Pt() < part2.Pt()) llep = lep2->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");

      histAddVal(cos(absnormPhi(part1.Phi() - _MET->phi())), "Part1CosDphiPtandMet");
      histAddVal(cos(absnormPhi(part2.Phi() - _MET->phi())), "Part2CosDphiPtandMet");


      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }

      double InvMass = diParticleMass(part1,part2, "InvariantMass");
      histAddVal(InvMass, "InvariantMass");

      double diMass1 = CalculateDiLepMassDeltaPt(part1, part2, mass1, mass2);
      histAddVal(diMass1, "ReconstructableMassDeltaPt");

      double cosDphiLeadLepMet = cos(absnormPhi(llep.Phi() - _MET->phi()));
      histAddVal(cosDphiLeadLepMet, "RecoCosDphiPtLeadLepandMet");

      double AbscosDphiLeadLepMet = abs(cos(absnormPhi(llep.Phi() - _MET->phi())));
      histAddVal(AbscosDphiLeadLepMet, "RecoAbsCosDphiPtLeadLepandMet");

  	  if(calculateLeptonMetMt(llep) > 0){
  	  	histAddVal(diMass1, "ReconstructableMassDeltaPt_nonzeroMt");
  	  	histAddVal(cosDphiLeadLepMet, "RecoCosDphiPtLeadLepandMet_nonzeroMt");
  	  	histAddVal(AbscosDphiLeadLepMet, "RecoAbsCosDphiPtLeadLepandMet_nonzeroMt");
  	  }

      double ptSum = part1.Pt() + part2.Pt();
      histAddVal(ptSum, "SumOfPt");

      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt");
      histAddVal(lep2->charge(p2) * lep1->charge(p1), "OSLS");
      histAddVal(PZeta, "PZeta");
      histAddVal(PZetaVis, "PZetaVis");
      histAddVal2(PZetaVis,PZeta, "Zeta2D");
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      // Diparticle pT
      histAddVal((part1 + part2).Pt(), "Pt");

      // DeltaPhi variables
      float minDPhi1Met = 9999.9, minDPhi2Met = 9999.9;

      if(active_part->at(CUTS::eRJet1)->size() > 0){

        for(size_t i = 0; i < active_part->at(CUTS::eRJet1)->size(); i++){
          TLorentzVector jet1P4 = _Jet->p4(active_part->at(CUTS::eRJet1)->at(i));
          float dPhi1Met = normPhi(jet1P4.Phi() - _MET->phi());
          float absdPhi1Met = absnormPhi(jet1P4.Phi() - _MET->phi());

          histAddVal(dPhi1Met, "Jet1MetDPhi");
          histAddVal(absdPhi1Met, "Jet1MetAbsDPhi");

          if(lep2->charge(p2) * lep1->charge(p1) < 0){
            histAddVal(dPhi1Met, "Jet1MetDPhiOS");
            histAddVal(absdPhi1Met, "Jet1MetAbsDPhiOS");
          }
          else if(lep2->charge(p2) * lep1->charge(p1) > 0){
            histAddVal(dPhi1Met, "Jet1MetDPhiLS");
            histAddVal(absdPhi1Met, "Jet1MetAbsDPhiLS");
          }

          if(absdPhi1Met < minDPhi1Met){
            minDPhi1Met = absdPhi1Met;
          }
        }

        histAddVal(minDPhi1Met, "Jet1MetMinAbsDPhi");
        if(lep2->charge(p2) * lep1->charge(p1) < 0){ histAddVal(minDPhi1Met, "Jet1MetMinAbsDPhiOS"); }
        else if(lep2->charge(p2) * lep1->charge(p1) > 0){ histAddVal(minDPhi1Met, "Jet1MetMinAbsDPhiLS"); }
      }

      if(active_part->at(CUTS::eRJet2)->size() > 0){

        for(size_t i = 0; i < active_part->at(CUTS::eRJet2)->size(); i++){
          TLorentzVector jet2P4 = _Jet->p4(active_part->at(CUTS::eRJet2)->at(i));
          float dPhi2Met = normPhi(jet2P4.Phi() - _MET->phi());
          float absdPhi2Met = absnormPhi(jet2P4.Phi() - _MET->phi());

          histAddVal(dPhi2Met, "Jet2MetDPhi");
          histAddVal(absdPhi2Met, "Jet2MetAbsDPhi");

          if(lep2->charge(p2) * lep1->charge(p1) < 0){
            histAddVal(dPhi2Met, "Jet2MetDPhiOS");
            histAddVal(absdPhi2Met, "Jet2MetAbsDPhiOS");
          }
          else if(lep2->charge(p2) * lep1->charge(p1) > 0){
            histAddVal(dPhi2Met, "Jet2MetDPhiLS");
            histAddVal(absdPhi2Met, "Jet2MetAbsDPhiLS");
          }

          if(absdPhi2Met < minDPhi2Met){
            minDPhi2Met = absdPhi2Met;
          }
        }

        histAddVal(minDPhi2Met, "Jet2MetMinAbsDPhi");
        if(lep2->charge(p2) * lep1->charge(p1) < 0){ histAddVal(minDPhi2Met, "Jet2MetMinAbsDPhiOS"); }
        else if(lep2->charge(p2) * lep1->charge(p1) > 0){ histAddVal(minDPhi2Met, "Jet2MetMinAbsDPhiLS"); }
      }

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }

      if(lep1->type != PType::Tau) {
        histAddVal(isZdecay(part1, *lep1), "Part1IsZdecay");
      }
      if(lep2->type != PType::Tau){
        histAddVal(isZdecay(part2, *lep2), "Part2IsZdecay");
      }
    }
  }
}

void Analyzer::initializePileupInfo(const bool& specialPU, std::string outfilename){
   // If specialPUcalculation is true, then take the name of the output file (when submitting jobs)
   // to retrieve the right MC nTruePU distribution to calculate the PU weights, otherwise, it will do
   // the calculation as usual, having a single MC nTruePU histogram (MCHistos option).
   // If you need to run the Analyzer interactively and use the specialPUcalculation option, make sure that you
   // include the name of the sample to analyze. For example, you can look at the names given to output files in past
   // runs sumbitted as jobs and use those names instead. This will be improved later.

   // try-catch block for initialize pileup info
   try{

     if(!specialPU){// No special PU calculation.
       initializePileupWeights(distats["Run"].smap.at("MCHistos"),distats["Run"].smap.at("DataHistos"),distats["Run"].smap.at("DataPUHistName"),distats["Run"].smap.at("MCPUHistName"));
     }
     else{ // Special PU calculation
    std::string outputname = outfilename;

    std::string delimitertune = "_Tune";
    std::string delimiterenergy = "_13TeV";

    bool istuneinname = outputname.find(delimitertune.c_str()) != std::string::npos;

    std::string samplename;

    if((samplename.length() == 0) && istuneinname){
      unsigned int pos = outputname.find(delimitertune.c_str());
      samplename = outputname.erase(pos, (outputname.substr(pos).length()));
    }
    else if((samplename.length() == 0) && (!istuneinname)){
      unsigned int pos = outputname.find(delimiterenergy.c_str());
      samplename = outputname.erase(pos, (outputname.substr(pos).length()));
    }

         initializePileupWeights(distats["Run"].smap.at("SpecialMCPUHistos"),distats["Run"].smap.at("DataHistos"),distats["Run"].smap.at("DataPUHistName"), samplename);
     }
   }// end of try block
   catch(std::runtime_error& err){
    std::cerr << "ERROR in initializePileupInfo! " << std::endl;
    std::cerr << "\t" << err.what()  << std::endl;
    std::cout << "\tMake sure you are using the right file and the pileup histogram name is correct in Run_info.in." << std::endl;
    if(specialPUcalculation){
      std::cerr << "\tWARNING: You are using SpecialPUCalculation. Check that the output filename (" << outfilename << ") contains the name of the analyzed MC sample." << std::endl;
    }
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   }
    catch(std::out_of_range& err){
    std::cerr << "ERROR in initializePileupInfo! " << std::endl;
    std::cerr << "\t" << err.what()  << std::endl;
    if(specialPU){
      std::cerr << "\tYou are using SpecialPUCalculation. Name of the MC sample was not found in the output filename (" << outfilename << ")." << std::endl;
      std::cerr << "\tTry with a different name that contains the name of your MC sample followed by _Tune or _13TeV (e.g. DYJetsToLL_M-50_TuneCP5.root)." << std::endl;
    }
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   }
   catch(...){
    std::cerr << "ERROR in initializePileupInfo! Unknown exception." << std::endl;
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   } // end of catch blocks

 }


void Analyzer::initializePileupWeights(std::string MCHisto, std::string DataHisto, std::string DataHistoName, std::string MCHistoName) {

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny(MCHistoName.c_str());
  if(!histmc) throw std::runtime_error(("Failed to extract histogram "+MCHistoName+" from "+PUSPACE+MCHisto+"!").c_str());
  histmc->Scale(1./histmc->Integral());
  histmc->SetDirectory(0);

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  hist_pu_wgt = dynamic_cast<TH1D*> (file2->FindObjectAny(DataHistoName.c_str()) );
  if(!hist_pu_wgt) throw std::runtime_error(("Failed to extract histogram "+DataHistoName+" from "+PUSPACE+DataHisto+"!").c_str());

  hist_pu_wgt->Scale(1.0/hist_pu_wgt->Integral());
  hist_pu_wgt->Divide(histmc);
  hist_pu_wgt->SetDirectory(0);

  hist_pu_wgt_up = dynamic_cast<TH1D*> ( file2->FindObjectAny((DataHistoName+"Up").c_str()) );
  if(hist_pu_wgt_up){
    hist_pu_wgt_up->Scale(1.0/hist_pu_wgt_up->Integral());
    hist_pu_wgt_up->Divide(histmc);
    hist_pu_wgt_up->SetDirectory(0);
  }

  hist_pu_wgt_do = dynamic_cast<TH1D*> ( file2->FindObjectAny((DataHistoName+"Do").c_str()) );
  if(hist_pu_wgt_do){
    hist_pu_wgt_do->Scale(1.0/hist_pu_wgt_do->Integral());
    hist_pu_wgt_do->Divide(histmc);
    hist_pu_wgt_do->SetDirectory(0);
  }
  /*
  //double factor = histmc->Integral() / histdata->Integral();
  float value,valueUp,valueDown;
  // The bin corresponding to nTruePU = 0 will be bin = 1 and, the calculated weight will be stored
  // in hPU[bin]. That's why when we calculate the pu_weight in setupEventGeneral, we require hPU[(int)nTruePU+1].
  for(int bin=0; bin <= histmc->GetNbinsX(); bin++) {
    if(histmc->GetBinContent(bin) == 0){
      value = 1;
      valueUp = 1;
      valueDown = 1;
    }else{
      value = histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
      std::cout << "bin = " << bin << ", value pu wgt = " << value << std::endl;
      if(histdata_up){
        valueUp = histdata_up->GetBinContent(bin) / histmc->GetBinContent(bin);
        valueDown = histdata_down->GetBinContent(bin) / histmc->GetBinContent(bin);
      }
    }
    hPU[bin]      = value;
    if(histdata_up){
      hPU_up[bin]   = valueUp;
      hPU_down[bin] = valueDown;
    }else{
      hPU_up[bin]   = value;
      hPU_down[bin] = value;
    }
  }
  std::cout << " ---- " << std::endl;
  */
  /*
  for(int bin=0; bin<= hist_pu_wgt->GetNbinsX(); bin++){
    std::cout << "bin = " << bin << ", pu wgt (hist) = " << ( hist_pu_wgt->GetBinContent(bin) ) << std::endl;
  }
  */
  file1->Close();
  file2->Close();

}

void Analyzer::initializeNPVWeights(std::string year){

  TFile *file = new TFile((PUSPACE+"NPVweights.root").c_str());
  histnpvwgt = dynamic_cast<TH1D*>( file->FindObjectAny( ("NVertices_"+year).c_str() ) );
  if(!histnpvwgt) throw std::runtime_error(("Failed to extract NPV weight histogram from "+PUSPACE+"NPVweights.root!").c_str());
  histnpvwgt->SetDirectory(0);
  // TH1D* histnpvwgt_up = (TH1D*)file1->FindObjectAny((DataHistoName+"Up").c_str());
  // TH1D* histnpvwgt_down = (TH1D*)file1->FindObjectAny((DataHistoName+"Do").c_str());
  //std::cout << "Initialized NPV weights successfully." << std::endl;
  file->Close();
}

void Analyzer::initializeWkfactor(std::vector<std::string> infiles) {
  if(infiles[0].find("WJets") != std::string::npos){
    isWSample = true;
  }else{
    isWSample=false;
    return;
  }
  //W-jet k-factor Histograms:
  TFile k_ele("Pileup/k_faktors_ele.root");
  TFile k_mu("Pileup/k_faktors_mu.root");
  TFile k_tau("Pileup/k_faktors_tau.root");

  k_ele_h =dynamic_cast<TH1D*>(k_ele.FindObjectAny("k_fac_m"));
  k_mu_h  =dynamic_cast<TH1D*>(k_mu.FindObjectAny("k_fac_m"));
  k_tau_h =dynamic_cast<TH1D*>(k_tau.FindObjectAny("k_fac_m"));

  k_ele.Close();
  k_mu.Close();
  k_tau.Close();

}

///Normalizes phi to be between -PI and PI
double normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}


///Takes the absolute value of of normPhi (made because constant use)
double absnormPhi(double phi) {
  return abs(normPhi(phi));
}

#include "BTaggingSFProd.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

void BTaggingSFProd::init(PartStats& stats, std::string pathtoeffmap, std::string year){

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

	// -------- b-tagging scale factors setup --------- //

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
	 btagcalib = BTagCalibration(btagalgoname, (pathtoeffmap+"BJetDatabase/"+btagsffilename).c_str());

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
	 std::cout << "\t Setting dummy b-tagging from " << (pathtoeffmap+"btagging.csv").c_str() << std::endl;
	 btagcalib = BTagCalibration("DeepCSV", (pathtoeffmap+"btagging.csv").c_str());
 }

// Check BTagCalibrationStandalone.h for more information

	static std::map<std::string, int> bjetflavors = {
	 {"bjet", 0}, {"cjet", 1}, {"lightjet", 2},
	};

	bjetflavor = (BTagEntry::JetFlavor) bjetflavors["bjet"];

	static std::map<std::string, int> btagoperatingpoints = {
	 {"loose", 0}, {"medium", 1}, {"tight", 2}, {"reshaping", 3}
	};

	std::string btagwp = stats.smap.at("JetBTaggingWP").c_str();

	b_workingpoint = (BTagEntry::OperatingPoint) btagoperatingpoints[btagwp];

	// ----------------- //

	// -------- b-tagging efficiency setup --------- //
	// Set the ROOT file which contains the b-tagging efficiency for your topology
	effcyMap_filename = Form("%sBTaggingEffcy.root", pathtoeffmap.data());

	// Set the histogram name from where the efficiency will be reader
	effcyMap_histoName = Form("h2_eff_mc_%s_%s_%s", btagalgoname.data(), btagwp.data(), year.data());

	// Create file for efficiency map
	TFile* effcyMap_file = new TFile(effcyMap_filename);
	if(!effcyMap_file || effcyMap_file->IsZombie()){
		std::cerr << "(BTaggingSFProd) ERROR! Failed to open input file " << effcyMap_filename << std::endl;
		std::cerr << "(BTaggingSFProd) Check that this file is in the specified path: " << pathtoeffmap << std::endl;
		std::cerr << "\tAborting Analyzer..." << std::endl;
		std::abort();
	}

	// Get the efficiency map histogram
	effcyMap_histo = dynamic_cast<TH2F*>(effcyMap_file->Get(effcyMap_histoName));
	effcyMap_histo->SetDirectory(0);
	if(!effcyMap_histo) throw std::runtime_error(("(BTaggingSFProd) Failed to extract histogram "+effcyMap_histoName+" from "+effcyMap_filename+"!"));
	effcyMap_file->Close();
	delete effcyMap_file;
	// ----------------- //

}

float BTaggingSFProd::calculateBTagSF(Jet& jets, GenJets& genjets, std::vector<int> passingbjets, std::vector<int> failingbjets, const std::string systname){

	// Load the info from the btaggin reader
  btagsfreader = BTagCalibrationReader(b_workingpoint, systname);
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  int matchedGenJetIndex = -1;
  int jetPartonFlavor = 0;
  int jetHadronFlavor = 0;

  for(size_t i = 0; i < passingbjets.size(); i++){

		matchedGenJetIndex = jets.genJetIdx[passingbjets[i]];
		jetPartonFlavor = abs(genjets.genPartonFlavor[matchedGenJetIndex]);
		jetHadronFlavor = static_cast<unsigned>(genjets.genHadronFlavor[matchedGenJetIndex]);

		// Check that this jet is a genuine b-jets
		if( (abs(jetPartonFlavor) != 5) || (abs(jetHadronFlavor) != 5) ) continue;

		// Get the Lorentz vector for the corresponding jet
    TLorentzVector bjetP4 = jets.p4(passingbjets[i]);

		// Get the SF and efficiency for this b-jet.
		float bjetSF = btagsfreader.eval_auto_bounds(systname, bjetflavor, bjetP4.Eta(), bjetP4.Pt());
		float bjetEff = getBJetEfficiency(bjetP4.Eta(), bjetP4.Pt(), effcyMap_histo);

		// Push the values for the corresponding probabilities
		float prob_mc = bjetEff;
		float prob_data = bjetSF*bjetEff;

		P_BJets_MC.push_back(prob_mc);
		P_BJets_Data.push_back(prob_data);
	}

  for(size_t i = 0; i < failingbjets.size(); i++){

		matchedGenJetIndex = jets.genJetIdx[passingbjets[i]];
		jetPartonFlavor = abs(genjets.genPartonFlavor[matchedGenJetIndex]);
		jetHadronFlavor = static_cast<unsigned>(genjets.genHadronFlavor[matchedGenJetIndex]);

		// Check that this jet is a genuine b-jets
		if( (abs(jetPartonFlavor) != 5) || (abs(jetHadronFlavor) != 5) ) continue;

		// Get the Lorentz vector for the corresponding jet
    TLorentzVector bjetP4 = jets.p4(failingbjets[i]);

		// Get the SF and efficiency for this b-jet.
		float bjetSF = btagsfreader.eval_auto_bounds(systname, bjetflavor, bjetP4.Eta(), bjetP4.Pt());
		float bjetEff = getBJetEfficiency(bjetP4.Eta(), bjetP4.Pt(), effcyMap_histo);

		// Push the values for the corresponding probabilities
		float prob_mc = 1.0 - bjetEff;
		float prob_data = 1.0 - bjetSF*bjetEff;

		P_BJets_MC.push_back(prob_mc);
		P_BJets_Data.push_back(prob_data);

	}

  float bjetSF = produceBJetSF(P_BJets_Data, P_BJets_MC);

  P_BJets_Data.clear();
	P_BJets_Data.shrink_to_fit();
  P_BJets_MC.clear();
	P_BJets_MC.shrink_to_fit();

  return bjetSF;
}

float BTaggingSFProd::getBJetEfficiency(double eta, double pt, TH2F* h_effmap){

	if(h_effmap == nullptr){
		std::cout << "B-tagging efficiency map not found, setting b-tagging efficiency to 1.0." << std::endl;
		return 1.0;
	}

	// Check pt is not above map overflow
	int nbinsy = h_effmap->GetNbinsY();
	float maxy = h_effmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
	if(pt >= maxy)
		pt = maxy - 0.01;
	int thebin = h_effmap->FindBin(eta,pt);

	float efficiency = h_effmap->GetBinContent(thebin);

	return efficiency;

}

float BTaggingSFProd::produceBJetSF(std::vector<float> Prob_Data, std::vector<float> Prob_MC){

	if(!(Prob_Data.size() == Prob_MC.size())){
	 std::cerr << "(BTaggingSFProd) ERROR! Probability b-tagging SF arrays have different sizes! Returning 1.0." << std::endl;
	 return 1.0;
	}

	float total_prob_Data = 1.0;
	float total_prob_MC = 1.0;

	for(size_t i = 0; i < Prob_Data.size(); i++){

		total_prob_Data *= Prob_Data[i];
		total_prob_MC *= Prob_MC[i];
	}

	if(total_prob_MC == 0){
		std::cerr << "(BTaggingSFProd) ERROR! Probability b-tagging SF for MC is zero, cannot divide by zero! Returning 1.0." << std::endl;
		return 1.0;
	}

	float btaggingSF = total_prob_Data/total_prob_MC;

	return btaggingSF;
}

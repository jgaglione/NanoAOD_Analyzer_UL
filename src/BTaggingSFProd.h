#ifndef BTaggingSFProd_h
#define BTaggingSFProd_h

/* -*- C++ -*-
 *
 * Class: L1ECALPrefiringWgtProd
 * Based on recipe by Laurent Thomas @ CMS
 * Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
 * Source: https://github.com/cms-sw/cmssw/blob/8706dbe8a09e7e1314f2127288cfc39051851eea/PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc
 * Goal: calculate event weights for 2016+2017 to mitigate the L1 ECAL prefiring issue.
 * Adapted by: Brenda Fabela.
 * Date: Sept 10, 2020
*/
#include <iostream>
#include <memory>
#include <TFile.h>   // TFile
#include <TString.h> // Form
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <TH2.h>
#include "Particle.h"
#include "./btagging/BTagCalibrationStandalone.h"

enum variations {central = 0, up, down};

class BTaggingSFProd{
	public:
		BTaggingSFProd() { };
		~BTaggingSFProd() { };
		void init(PartStats& stats, std::string pathtoeffmap, std::string year);
		float calculateBTagSF(Jet& jets, GenJets& genjets, std::vector<int> passingbjets, std::vector<int> failingbjets, const std::string systname);

	private:

		BTagEntry::JetFlavor bjetflavor;
	  BTagEntry::OperatingPoint b_workingpoint;

		// B-tagging scale factors - calibration + readers
	  BTagCalibration btagcalib;
	  BTagCalibrationReader btagsfreader;

		std::vector<float> P_BJets_Data;
		std::vector<float> P_BJets_MC;

		float getBJetEfficiency(double eta, double pt, TH2F* h_effmap);
		float produceBJetSF(std::vector<float> Prob_Data, std::vector<float> Prob_MC);

		TString effcyMap_filename;
		TString effcyMap_histoName;

		TH2F* effcyMap_histo;

};

#endif

#include <cstring>
#include <iostream>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include "Tracking.h"

// void PrintHelp() {
// 	std::cout << "usage: ./Tracking [options]" << std::endl;
// 	std::cout << "options:" << std::endl;
// 	std::cout << "  -h 		print this information" << std:: endl;
// 	std::cout << "  -t 		track, create tracking.root" << std::endl;
// 	std::cout << "  -e 		effciency, create effciency.root" << std::endl;
// 	std::cout << "  -a 		track and effciency, -e and -t" << std::endl;
// 	exit(0);
// }

// bool doTrack = false;
// bool doEffciency = false;


int main(int argc, char **argv) {
	// if (argc != 2 || argv[1][0] != '-' || strlen(argv) != 2) {
	// 	PrintHelp();
	// }
	// switch (argv[1][1]) {
	// 	case 't':
	// 		doTrack = true;
	// 		break;
	// 	case 'e':
	// 		doEffciency = true;
	// 		break;
	// 	case 'a':
	// 		doTrack = true;
	// 		doEffciency = true;
	// 		break;
	// 	case 'h':
	// 	default:
	// 		PrintHelp();
	// }

	if (argc != 1) {
		std::cout << "usage: ./Tracking" << std::endl;
	}

	TString inFile("../data/f8ppac001.root");
	TFile *ipf = new TFile(inFile.Data());
	if (!ipf->IsOpen()) {
		std::cout << "Error open file " << inFile << std::endl;
	}
	TTree *ipt = (TTree*)ipf->Get("tree");

	Tracking *tk = new Tracking(ipt);
	tk->Loop();

	ipf->Close();
	return 0;
}
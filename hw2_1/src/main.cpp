#include <cstring>
#include <iostream>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <ctime>
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

	clock_t t;

	TString inFile("../data/f8ppac001.root");
	TFile *ipf = new TFile(inFile.Data());
	if (!ipf->IsOpen()) {
		std::cout << "Error open file " << inFile << std::endl;
	}
	TTree *ipt = (TTree*)ipf->Get("tree");

	TFile *opf;
	TTree *opt;

	// opf = new TFile("../data/tracking.root", "recreate");
	// opt = new TTree("tree", "ppac tracking traditional fit method");
	// Tracking *tk = new Tracking(ipt);
	// t = clock();
	// tk->Loop(opt, false);
	// t = clock() - t;
	// printf("Old loop time: %f s\n", float(t) / CLOCKS_PER_SEC);
	// opt->Write();
	// opf->Close();

	opf = new TFile("../data/tracking.root", "recreate");
	opt = new TTree("tree", "ppac tracking new fit method");
	Tracking *tkf = new Tracking(ipt);
	t = clock();
	tkf->Loop(opt);
	t = clock() - t;
	printf("New loop time: %f s\n", float(t) / CLOCKS_PER_SEC);
	opt->Write();
	opf->Close();
	// delete opt 		no need to do this, tree was deleted when file closed



	ipf->Close();

	return 0;
}
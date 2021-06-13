#include <TFile.h>

#include "DSSD.h"


int main(int argc, char **argv) {
	TFile *ipf = new TFile("../data/data_16Cs.root", "read");
	TTree *tree = (TTree*)ipf->Get("tree");

	TFile *opf = new TFile("../data/norm_16C.root", "recreate");
	TTree *opt = new TTree("tree", "noramlize tree for 16C");

	TFile *logFile = new TFile("../data/normLog.root", "recreate");

	TFile *graphFile = new TFile("../data/graph_16C.root", "update");

	MDSSD *md = new MDSSD(tree, opt);
	md->SetLogFile(logFile);
	md->SetGraphFile(graphFile);

	md->AddDSSD("d1", 32, 32);
	md->AddDSSD("d2", 32, 32);
	md->AddDSSD("d3", 32, 32);

	md->AddNode(0, 16, 17, 9, 16, DIR_X2Y);
	md->AddNode(0, 11, 24, 9, 16, DIR_Y2X);
	md->AddNode(0, 11, 24, 0, 32, DIR_X2Y);
	md->AddNode(0, 0, 32, 0, 32, DIR_Y2X);

	md->AddNode(1, 16, 17, 8, 18, DIR_X2Y);
	md->AddNode(1, 10, 24, 8, 18, DIR_Y2X);
	md->AddNode(1, 10, 24, 0, 32, DIR_X2Y);
	md->AddNode(1, 0, 32, 0, 32, DIR_Y2X);

	md->AddNode(2, 15, 16, 8, 18, DIR_X2Y);
	md->AddNode(2, 8, 24, 8, 18, DIR_Y2X);
	md->AddNode(2, 8, 24, 0, 32, DIR_X2Y);
	md->AddNode(2, 0, 32, 0, 32, DIR_Y2X);

	md->Norm();

	opf->cd();
	opt->Write("tree", TObject::kOverwrite);
	opf->Close();
	ipf->Close();
	logFile->Close();
	graphFile->Close();

	return 0;
}
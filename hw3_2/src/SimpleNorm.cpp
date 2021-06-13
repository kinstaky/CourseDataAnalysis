#include <TFile.h>

#include "DSSD.h"


int main(int argc, char **argv) {
	TFile *ipf = new TFile("../data/data_16Cs.root", "read");
	TTree *tree = (TTree*)ipf->Get("tree");

	TFile *opf = new TFile("../data/normS_16C.root", "recreate");
	TTree *opt = new TTree("tree", "noramlize tree for 16C");

	TFile *logFile = new TFile("../data/normSLog.root", "recreate");

	TFile *graphFile = new TFile("../data/graph_16C.root", "update");

	MDSSD *md = new MDSSD(tree, opt);
	md->SetLogFile(logFile);
	md->SetGraphFile(graphFile);

	md->AddDSSD("d1", 32, 32);
	md->AddDSSD("d2", 32, 32);
	md->AddDSSD("d3", 32, 32);

	md->AddNode(0, 15, 16, 0, 32, DIR_X2Y);
	md->AddNode(0, 0, 32, 12, 13, DIR_Y2X);
	md->AddNode(1, 15, 16, 0, 32, DIR_X2Y);
	md->AddNode(1, 0, 32, 15, 16, DIR_Y2X);
	md->AddNode(2, 15, 16, 0, 32, DIR_X2Y);
	md->AddNode(2, 0, 32, 15, 16, DIR_Y2X);

	md->Norm();

	opf->cd();
	opt->Write("tree", TObject::kOverwrite);
	opf->Close();
	ipf->Close();
	logFile->Close();

	return 0;
}
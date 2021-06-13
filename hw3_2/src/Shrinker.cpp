#include <TFile.h>
#include <TTree.h>

#include "Shrink.h"

int main(int argc, char **argv) {
	TFile *ipf = new TFile("../data/norm_16C.root", "read");
	TTree *ipt = (TTree*)ipf->Get("tree");

	TFile *sf = new TFile("../data/shrink_16C.root", "recreate");
	TTree *st = new TTree("tree", "shrink tree of 16C");

	TFile *rf = new TFile("../data/rest_16C.root", "recreate");
	TTree *rt = new TTree("tree", "rest tree for shrinker 16C");

	DSSDShrink *ds = new DSSDShrink(ipt, st, rt);

	ds->AddDSSD("d1", 32, 32, 70);
	ds->AddDSSD("d2", 32, 32, 50);
	ds->AddDSSD("d3", 32, 32, 50);

	ds->Shrink();

	sf->cd();
	st->Write();
	sf->Close();

	rf->cd();
	rt->Write();
	rf->Close();

	ipf->Close();
	return 0;
}
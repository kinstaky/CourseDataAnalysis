#include "dssd.h"

#include <TTree.h>
#include <TFile.h>

int main() {
	TFile *ipf = new TFile("../data/alpha.root", "read");
	TTree *ipt = (TTree*)ipf->Get("tree");

	TFile *opf = new TFile("../data/alphaC.root", "recreate");
	TTree *opt = new TTree("tree", "dssd after front-back correlation");

	DSSD* d = new DSSD(ipt, opt, 48, 128);

	d->Analysis();

	opf->Close();
	ipf->Close();
	return 0;
}
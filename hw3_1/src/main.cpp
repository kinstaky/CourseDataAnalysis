#include "Adssd.h"

const int alphaNum = 4;
const Double_t alphaPeaks[alphaNum] {5.5658, 6.1748, 6.6708, 8.6931};

int main(int argc, char **argv) {
	TFile *ipf = new TFile("../data/s4c.root", "read");
	TTree *ipt = (TTree*)ipf->Get("tree");

	TFile *opf = new TFile("../data/s4calibration.root", "recreate");
	TTree *opt = new TTree("tree", "alpha calibration");

	TFile *logf = new TFile("../data/s4Log.root", "recreate");

	Adssd *ad = new Adssd(ipt, logf);
	ad->SetAlpha(alphaNum, alphaPeaks);

	ad->Analysis(opt);

	opf->cd();
	opt->Write();
	opf->Close();
	logf->Close();
	ipf->Close();

	return 0;
}
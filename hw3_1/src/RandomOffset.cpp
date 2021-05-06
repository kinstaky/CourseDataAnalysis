#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TRandom3.h>


const int chp = 48;
const int chr = 48;


int main(int argc, char **argv) {
	TFile *ipf = new TFile("../data/s4.root");
	TTree *ipt = (TTree*)ipf->Get("tree");

	Int_t pe[chp], re[chr];

	ipt->SetBranchAddress("pe", pe);
	ipt->SetBranchAddress("re", re);


	TFile *opf = new TFile("../data/s4c.root", "recreate");
	TTree *opt = new TTree("tree", "alpha with random offset");

	TRandom3 *gr = new TRandom3(0);
	Double_t cpe[chp], cre[chr];
	opt->Branch("pe", pe, TString::Format("pe[%d]/I", chp));
	opt->Branch("re", re, TString::Format("re[%d]/I", chr));
	opt->Branch("cpe", cpe, TString::Format("cpe[%d]/D", chp));
	opt->Branch("cre", cre, TString::Format("cre[%d]/D", chr));

	Long64_t nentries = ipt->GetEntries();
	Long64_t nentries100 = nentries / 100;
	printf("randomOffset   0%%");
	for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
		ipt->GetEntry(jentry);
		for (int i = 0; i != chp; ++i) {
			cpe[i] = pe[i] + gr->Uniform(0, 1);
		}
		for (int i = 0; i != chr; ++i) {
			cre[i] = re[i] + gr->Uniform(0, 1);
		}
		opt->Fill();

		if ((jentry+1) % (nentries100) == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1) / nentries100);
			fflush(stdout);
		}
	}
	printf("\n");

	opt->Write();
	opf->Close();
	ipf->Close();
	return 0;
}
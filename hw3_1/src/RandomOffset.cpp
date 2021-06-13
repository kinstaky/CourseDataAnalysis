#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TString.h>
#include <TRandom3.h>
#include <TH1F.h>


const int chp = 48;
const int chr = 48;

TH1F *hp[chp];
TH1F *hr[chr];
TH1F *hpAll, *hrAll;

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

	hpAll = new TH1F("hpe", "hpe", 2000, 0, 2000);
	hrAll = new TH1F("hre", "hre", 2000, 0, 2000);
	TDirectoryFile *dph = new TDirectoryFile("phist", "phist");
	dph->cd();
	for (int i = 0; i != chp; ++i) {
		TString hName;
		hName.Form("hpe%d", i);
		hp[i] = new TH1F(hName.Data(), hName.Data(), 2000, 0, 2000);
	}
	opf->cd();
	TDirectoryFile *drh = new TDirectoryFile("rhist", "rhist");
	for (int i = 0; i != chr; ++i) {
		TString hName;
		hName.Form("hre%d", i);
		hr[i] = new TH1F(hName.Data(), hName.Data(), 2000, 0, 2000);
	}

	Long64_t nentries = ipt->GetEntries();
	Long64_t nentries100 = nentries / 100;
	printf("RandomOffset   0%%");
	for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
		ipt->GetEntry(jentry);
		for (int i = 0; i != chp; ++i) {
			cpe[i] = pe[i] + gr->Uniform(0, 1);
			hp[i]->Fill(cpe[i]);
			hpAll->Fill(cpe[i]);
		}
		for (int i = 0; i != chr; ++i) {
			cre[i] = re[i] + gr->Uniform(0, 1);
			hr[i]->Fill(cre[i]);
			hrAll->Fill(cre[i]);
		}
		opt->Fill();

		if ((jentry+1) % (nentries100) == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1) / nentries100);
			fflush(stdout);
		}
	}
	printf("\n");

	opf->cd();
	hpAll->Write();
	hrAll->Write();
	opt->Write("tree", TObject::kOverwrite);
	dph->cd();
	for (int i = 0; i != chp; ++i) {
		hp[i]->Write();
	}
	drh->cd();
	for (int i = 0; i != chr; ++i) {
		hr[i]->Write();
	}
	opf->Close();
	ipf->Close();
	return 0;
}
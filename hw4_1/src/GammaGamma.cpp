#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

int main() {
	// input data
	Int_t ahit;
	Int_t aid[16];
	Double_t ae[16], at[16];
	// input file and tree
	TFile *ipf = new TFile("../data/eurica.root", "read");
	TTree *ipt = (TTree*)ipf->Get("tree");
	ipt->SetBranchAddress("ahit", &ahit);
	ipt->SetBranchAddress("aid", aid);
	ipt->SetBranchAddress("ae", ae);
	ipt->SetBranchAddress("at", at);

	// output data
	Int_t achit;
	Int_t aidx[256], aidy[256];
	Double_t aex[256], aey[256], atx[256], aty[256];
	// oputput file and tree
	TFile *opf = new TFile("../data/ggc.root", "recreate");
	TTree *opt = new TTree("tree", "gamma-gamma coincidence tree");
	opt->Branch("achit", &achit, "achit/I");
	opt->Branch("aidx", aidx, "aidx[achit]/I");
	opt->Branch("aidy", aidy, "aidy[achit]/I");
	opt->Branch("aex", aex, "aex[achit]/D");
	opt->Branch("aey", aey, "aey[achit]/D");
	opt->Branch("atx", atx, "atx[achit]/D");
	opt->Branch("aty", aty, "aty[achit]/D");
	// output matrix
	TH2D *hgg = new TH2D("hgg", "gamma-gamma matrix", 1500, 0, 1500, 1500, 0, 1500);


	printf("processing   0%%");
	fflush(stdout);
	Long64_t nentry = ipt->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		ipt->GetEntry(jentry);
		achit = 0;
		for (int i = 0; i != ahit; ++i) {
			for (int j = 0; j != ahit; ++j) {
				if (i==j) continue;
				atx[achit] = at[i] + 2000.0 / sqrt(ae[i]) - 80.0;
				aty[achit] = at[j] + 2000.0 / sqrt(ae[j]) - 80.0;
				if (abs(atx[achit]-aty[achit]) > 1000) continue;
				aidx[achit] = aid[i];
				aidy[achit] = aid[j];
				aex[achit] = ae[i];
				aey[achit] = ae[j];

				if (abs(atx[achit]-aty[achit]) < 200) {
					hgg->Fill(aex[achit], aey[achit]);
				}
				achit++;
			}
		}
		if (achit) opt->Fill();

		if (jentry % nentry100) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}

	opt->Write();
	hgg->Write();
	opf->Close();
	ipf->Close();

	printf("\b\b\b\b100%%\n");
	return 0;
}
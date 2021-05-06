#include <TFile.h>
#include <TTree.h>

const int chp = 48;

int main(int argc, char **argv) {
	Double_t pe[chp];
	Double_t cpe[chp];
	Int_t phit, pid[chp];
	Double_t per[chp], cper[chp];


	TFile *ipf = new TFile("../data/s4c.root", "read");
	TTree *ipt = (TTree*)ipf->Get("tree");
	ipt->SetBranchAddress("cpe", pe);

	TFile *ipfc = new TFile("../data/s4calibration.root", "read");
	TTree *iptc = (TTree*)ipfc->Get("tree");
	iptc->SetBranchAddress("ccpe", cpe);

	TFile *opf = new TFile("../data/s4Hit.root", "recreate");
	TTree *opt = new TTree("tree", "alpha hit tree");
	opt->Branch("phit", &phit, "phit/I");
	opt->Branch("pid", pid, "pid[phit]/I");
	opt->Branch("pe", per, "pe[phit]/D");
	opt->Branch("cpe", cper, "cpe[phit]/D");

	Long64_t nentry = ipt->GetEntries();
	Long64_t nentry100 = nentry / 100;
	printf("filling   0%%");
	fflush(stdout);
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", jentry / nentry100);
			fflush(stdout);
		}
		ipt->GetEntry(jentry);
		iptc->GetEntry(jentry);
		phit = 0;
		for (int i = 0; i != chp; ++i) {
			if (pe[i] > 0) {
				pid[phit] = i;
				per[phit] = pe[i];
				cper[phit] = cpe[i];
				++phit;
			}
		}
		opt->Fill();
	}
	printf("\n");

	opt->Write();
	opf->Close();
	ipf->Close();
	ipfc->Close();
}
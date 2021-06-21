#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>

#include <map>


using namespace std;


// input data
Int_t hit;
Double_t ae[16];
ULong64_t at[16];
Double_t afStrip[16], abStrip[16];
Double_t afe[16], abe[16];
Double_t me[16];

struct pInfo {					// particle info
	Double_t Energy;
	Double_t FStrip, BStrip;
};

// maps
multimap<ULong64_t, pInfo> mapHeavy, mapDecay;

int setamBranch(TTree *tree) {
	tree->SetBranchAddress("hit", &hit);
	tree->SetBranchAddress("ae", ae);
	tree->SetBranchAddress("time", at);
	tree->SetBranchAddress("afStrip", afStrip);
	tree->SetBranchAddress("abStrip", abStrip);
	tree->SetBranchAddress("afe", afe);
	tree->SetBranchAddress("abe", abe);
	tree->SetBranchAddress("me", me);
	return 0;
}


int fillMap(TTree *tree) {
	mapHeavy.clear();
	mapDecay.clear();

	pInfo pi;
	printf("filling maps   0%%");
	fflush(stdout);
	Long64_t nentry = tree->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		tree->GetEntry(jentry);
		if (hit == 2) continue;
		pi.Energy = ae[0];
		pi.FStrip = afStrip[0];
		pi.BStrip = abStrip[0];
		if (me[0] < 0) {
			mapDecay.insert(make_pair(at[0], pi));
		} else {
			mapHeavy.insert(make_pair(at[0], pi));
		}

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");
	return 0;
}


int printMap(const multimap<ULong64_t, pInfo> &map, Long64_t entry = 20) {
	printf("%5s%15s%15s%15s%15s\n", "i", "time", "energy", "fStrip", "bStrip");
	Long64_t jentry = 0;
	for (auto it = map.begin(); it != map.end() && jentry != entry; ++it) {
		printf("%5lld%15lld%15lf%15lf%15lf\n", jentry, it->first, it->second.Energy, it->second.FStrip, it->second.BStrip);
		++jentry;
	}
	return 0;
}

int energyTime(ULong64_t tw) {
	TH2F *hdet = new TH2F("hdet", "decayEnergy vs decayTime", 400, Double_t(-tw), Double_t(tw), 300, 5e3, 8e3);
	for (auto ih = mapHeavy.begin(); ih != mapHeavy.end(); ++ih) {
		for (auto id = mapDecay.lower_bound(ih->first-tw); id != mapDecay.end(); ++id) {
			if (id->first >= ih->first+tw) break;
			Double_t dfs = abs(id->second.FStrip - ih->second.FStrip);			// delta front strip
			Double_t dbs = abs(id->second.BStrip - ih->second.BStrip);			// delta back strip
			if (dfs >= 1.0 || dbs >= 1.0) continue;
			Long64_t dt = id->first - ih->first;								// delta time
			if (id->second.Energy >= 5000 && id->second.Energy <= 8000) hdet->Fill(dt, id->second.Energy);
			// hdet->Fill(dt, id->second.Energy);
		}
	}

	hdet->Write();
	return 0;
}


// get decay time, with low energy and up energy to open the window
int decayTime(Double_t elow, Double_t eup, Double_t &dtime, Double_t &sdtime) {
	TH2F *hdet = (TH2F*)gDirectory->Get("hdet");
	hdet->GetYaxis()->SetRangeUser(elow, eup);
	TH1F *decayRa = (TH1F*)hdet->ProjectionX("decayRa");

	TF1 *fbkg = new TF1("fbkg", "pol0", -4e10, -1e10);
	decayRa->Fit(fbkg, "NQR+");
	Double_t p0 = fbkg->GetParameter(0);

	TF1 *fdecay = new TF1("fdecay", "[0]+[1]*exp(-x*log(2.0)/[2])", 1, 4e10);
	fdecay->FixParameter(0, p0);
	fdecay->SetParameter(2, 3e9);
	decayRa->Fit(fdecay, "NQR+");
	dtime = fdecay->GetParameter(2);
	sdtime = fdecay->GetParError(2);
	decayRa->GetListOfFunctions()->Add(fdecay);
	decayRa->GetListOfFunctions()->Add(fbkg);

	decayRa->Write();
	return 0;
}

int main() {
	TFile *amf = new TFile("../data/am.root", "read");
	TTree *amtree = (TTree*)amf->Get("tree");

	TFile *df = new TFile("../data/decay.root", "recreate");

	setamBranch(amtree);

	fillMap(amtree);

	// printMap(mapHeavy);
	// printMap(mapDecay);

	energyTime(4e10);					// generating decayEnergy vs decayTime distribution

	Double_t dtime, sdtime;
	decayTime(6960, 7080, dtime, sdtime);
	printf("half-life time of 210Ra is %lf +- %lfs\n", dtime/1e9, sdtime/1e9);

	df->Close();
	amf->Close();
	return 0;
}
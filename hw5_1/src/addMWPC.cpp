#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>

#include <map>

using namespace std;

// input data of alphaC.root
Int_t ahit;							// dssd hit
Double_t ae[16];					// dssd max energy
ULong64_t at[16];					// dssd time
Double_t afStrip[16], abStrip[16];	// dssd strips
Double_t afe[16], abe[16];			// dssd energy


// input data of mwpc.root
Double_t rme;						// mwpc energy
ULong64_t rmt;						// mwpc time


// output data
Int_t chit;							// hit
Double_t ce[16];					// dssd max energy
ULong64_t ct[16];					// dssd time
Double_t cfStrip[16], cbStrip[16];	// dssd strips
Double_t cfe[16], cbe[16];			// dssd energy
Double_t cme[16];					// mwpc energy


multimap<ULong64_t, Double_t> mapMWPC;





int setOptBranch(TTree *opt) {
	opt->Branch("hit", &chit, "hit/I");
	opt->Branch("ae", ce, "ae[hit]/D");
	opt->Branch("time", ct, "time[hit]/l");
	opt->Branch("afStrip", cfStrip, "afStrip[hit]/D");
	opt->Branch("abStrip", cbStrip, "abStrip[hit]/D");
	opt->Branch("afe", cfe, "afe[hit]/D");
	opt->Branch("abe", abe, "abe[hit]/D");
	opt->Branch("me", cme, "me[hit]/D");
	return 0;
}

int setAlphaCBranch(TTree *atree) {
	atree->SetBranchAddress("hit", &ahit);
	atree->SetBranchAddress("e", ae);
	atree->SetBranchAddress("time", at);
	atree->SetBranchAddress("fStrip", afStrip);
	atree->SetBranchAddress("bStrip", abStrip);
	atree->SetBranchAddress("fe", afe);
	atree->SetBranchAddress("be", abe);
	return 0;
}

int setMWPCBranch(TTree *mtree) {
	mtree->SetBranchAddress("energy", &rme);
	mtree->SetBranchAddress("timestamp", &rmt);
	return 0;
}


int fillMapMWPC(TTree *mtree) {
	mapMWPC.clear();
	printf("filling mwpc map   0%%");
	fflush(stdout);
	Long64_t nentry = mtree->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		mtree->GetEntry(jentry);
		mapMWPC.insert(make_pair(rmt, rme));

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");
	return 0;
}


int printMapMWPC(Long64_t entry = 20) {
	printf("mwpc map size %ld\n", mapMWPC.size());
	auto mi = mapMWPC.begin();
	printf("%10s%15s%15s\n", "i", "energy", "time");
	for (Long64_t jentry = 0; jentry != entry; ++jentry) {
		if (mi == mapMWPC.end()) break;
		printf("%10lld%15lf%15lld\n", jentry, mi->second, mi->first);
		++mi;
	}
	return 0;
}

int searchWindow(TTree *atree, ULong64_t tw = 100000) {
	TH1I *hdt = new TH1I("hdt", "alpha-mwpc distribution(ns)", 20000, -100000, 100000);
	printf("generating hdt   0%%");
	fflush(stdout);
	Long64_t nentry = atree->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		atree->GetEntry(jentry);

		for (int i = 0; i != ahit; ++i) {
			auto mi = mapMWPC.lower_bound(at[i]-tw);
			for (; mi != mapMWPC.end(); ++mi) {
				if (mi->first >= at[i]+tw) break;
				int dt = at[i] - mi->first;
				hdt->Fill(dt);
			}
		}

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");

	hdt->Write();
	return 0;
}

int awCorrelation(TTree *atree, TTree *opt, ULong64_t tw, ULong64_t toff) {
	TH1I *hemc = new TH1I("hemc", "energy of conincident mwpc", 3500, 0, 35000);
	TH1I *hfme = new TH1I("hfme", "found mwpc energy numbers", 10, 0, 10);
	printf("dssd mwpc correlation   0%%");
	fflush(stdout);
	Long64_t nentry = atree->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		atree->GetEntry(jentry);
		chit = ahit;
		for (int i = 0; i != chit; ++i) {
			ce[i] = ae[i];
			ct[i] = at[i];
			cfStrip[i] = afStrip[i];
			cbStrip[i] = abStrip[i];
			cfe[i] = afe[i];
			cbe[i] = abe[i];
			cme[i] = -1.0;
		}

		int foundme = 0;
		auto mi = mapMWPC.lower_bound(at[0]+toff-tw);
		for (; mi != mapMWPC.end(); ++mi) {
			if (mi->first >= at[0]+toff+tw) break;
			++foundme;
			for (int i = 0; i != chit; ++i) {
				cme[i] = mi->second;
			}
			opt->Fill();
			hemc->Fill(mi->second);
		}
		hfme->Fill(foundme);

		if (!foundme) {
			opt->Fill();
		}

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");

	hemc->Write();
	hfme->Write();
	opt->Write();
	return 0;
}


int main() {
	TFile *af = new TFile("../data/alphaC.root", "read");
	TTree *atree = (TTree*)af->Get("tree");
	setAlphaCBranch(atree);

	TFile *mf = new TFile("../data/mwpc.root", "read");
	TTree *mtree = (TTree*)mf->Get("tree");
	setMWPCBranch(mtree);

	TFile *opf = new TFile("../data/am.root", "recreate");
	TTree *opt = new TTree("tree", "tree of dssd and mwpc");
	setOptBranch(opt);

	fillMapMWPC(mtree);
	printMapMWPC();

	searchWindow(atree);

	awCorrelation(atree, opt, 100, 1200);

	opf->Close();
	mf->Close();
	af->Close();
	return 0;
}
#include "dssd.h"

#include <TH2D.h>

using std::make_pair;
using std::multimap;

DSSD::DSSD(TTree *iipt, TTree *oopt, Int_t fs, Int_t bs) {
	ipt = iipt;
	opt = oopt;
	fStrip = fs;
	bStrip = bs;
}

DSSD::~DSSD() {
}

int DSSD::Analysis() {
	setBranch();
	fillMap();
	// printMap();
	searchWindow();
	hecorr(100, -550);			// energy correlation
	getHit(200, 0);
	correlation(100, -550);
	return 0;
}


int DSSD::fillMap() {
	dssdInfo ds;
	printf("filling map   0%%");
	fflush(stdout);
	Long64_t nentry = ipt->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		ipt->GetEntry(jentry);
		ds.Energy = re;
		ds.Strip = rs;
		if (side == 0) mapfs.insert(make_pair(rt, ds));
		else if (side == 1) mapbs.insert(make_pair(rt, ds));

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");


	auto ift = mapfs.begin();
	auto ibt = mapbs.begin();
	TH2D *heRaw = new TH2D("heRaw", "front-back correlation raw", 3000, 0, 30000, 3000, 0, 30000);
	int msize = mapfs.size();
	msize = msize < (int)mapbs.size() ? msize : (int)mapbs.size();
	for (int i = 0; i != msize; ++i) {
		heRaw->Fill(ift->second.Energy, ibt->second.Energy);
		++ift;
		++ibt;
	}
	heRaw->Write();
	return 0;
}

int DSSD::printMap(int entry) {
	auto ift = mapfs.begin();
	auto ibt = mapbs.begin();
	printf("%4s%15s%15s%10s%10s%10s%10s\n", "i", "front-t", "back-t", "front-s", "back-s", "front-e", "back-e");
	for (int i = 0; i != entry; ++i) {
		printf("%4d%15llu%15llu%10d%10d%10.lf%10.lf\n", i, ift->first, ibt->first, ift->second.Strip, ibt->second.Strip, ift->second.Energy, ibt->second.Energy);
		++ift;
		++ibt;
	}
	return 0;
}


int DSSD::setBranch() {
	ipt->SetBranchAddress("energy", &re);
	ipt->SetBranchAddress("strip", &rs);
	ipt->SetBranchAddress("side", &side);
	ipt->SetBranchAddress("timestamp", &rt);


	opt->Branch("hit", &chit, "hit/I");
	opt->Branch("e", cEnergy, "e[hit]/D");
	opt->Branch("time", cTime, "time[hit]/l");
	opt->Branch("fStrip", cfStrip, "fStrip[hit]/D");
	opt->Branch("bStrip", cbStrip, "bStrip[hit]/D");
	opt->Branch("fe", cfe, "fe[hit]/D");
	opt->Branch("be", cbe, "be[hit]/D");
	return 0;
}


int DSSD::searchWindow(ULong64_t tw, ULong64_t toff) {
	TH1I *hdt = new TH1I("hdt", "ft-bt distribution(ns)", 20000, -100000, 100000);
	for (auto ift = mapfs.begin(); ift != mapfs.end(); ++ift) {
		auto ibt = mapbs.lower_bound(ift->first-toff-tw);
		while (ibt != mapbs.end()) {
			if (ibt->first >= ift->first-toff+tw) break;
			int dt = (ift->first-toff)-ibt->first;
			hdt->Fill(dt);
			ibt++;
		}
	}
	hdt->Write();
	return 0;
}


int DSSD::hecorr(int tw, int toff) {
	TH2D *hec = new TH2D("heCorr", "front-back energy correlation", 3000, 0, 30000, 3000, 0, 30000);
	for (auto ift = mapfs.begin(); ift != mapfs.end(); ++ift) {
		auto ibt = mapbs.lower_bound(ift->first-toff-tw);
		for (; ibt != mapbs.end(); ++ibt) {
			if (ibt->first >= ift->first-toff+tw) break;
			hec->Fill(ift->second.Energy, ibt->second.Energy);
		}
	}
	hec->Write();
	return 0;
}

int DSSD::getHit(int tw, int toff) {
	dssdInfo *ds[10];
	TH1I *hfhit = new TH1I("hfhit", "hit", 10, 0, 10);
	TH1I *hbhit = new TH1I("hbhit", "hit", 10, 0, 10);
	TH1I *hfdt = new TH1I("hfdt", "dt", 20, -200, 200);
	TH1I *hbdt = new TH1I("hbdt", "dt", 20, -200, 200);
	for (auto ift = mapfs.begin(); ift != mapfs.end(); ++ift) {
		int hit = 0;
		Double_t t0 = 0.0;
		auto jft = mapfs.lower_bound(ift->first-toff);
		for (; jft != mapfs.end(); ++jft) {
			if (jft->first >= ift->first-toff+tw) break;
			if (hit == 0) t0 = jft->first;
			else if (hit == 1) t0 = jft->first - t0;
			if (hit < 10) ds[hit] = &(jft->second);
			++hit;
		}
		if (hit == 2) hfdt->Fill(t0);
		hfhit->Fill(hit);
		for (int i = 0; i != hit && i != 10; ++i) {
			ds[i]->Hit = hit;
			ds[i]->Index = i;
		}
	}

	for (auto ibt = mapbs.begin(); ibt != mapbs.end(); ++ibt) {
		int hit = 0;
		Double_t t0 = 0.0;
		auto jbt = mapbs.lower_bound(ibt->first-toff);
		for (; jbt != mapbs.end(); ++jbt) {
			if (jbt->first >= ibt->first-toff+tw) break;
			if (hit == 0) t0 = jbt->first;
			else if (hit == 1) t0 = jbt->first - t0;
			if (hit < 10) ds[hit] = &(jbt->second);
			++hit;
		}
		if (hit == 2) hbdt->Fill(t0);
		hbhit->Fill(hit);
		for (int i = 0; i != hit && i != 10; ++i) {
			ds[i]->Hit = hit;
			ds[i]->Index = i;
		}
	}

	hfhit->Write();
	hbhit->Write();
	hfdt->Write();
	hbdt->Write();
	return 0;
}


int DSSD::correlation(int tw, int toff) {
	const int eThre = 50;
	multimap<Double_t, dsfbInfo> mapfbs;
	dsfbInfo dsfb;
	TH2D *hxyse = new TH2D("hxyse", "fse-bse : fse", 2000, -1000, 1000, 3000, 0, 30000);
	TH2D *hxyec = new TH2D("hxyec", "front-back correclation energy", 3000, 0, 30000, 3000, 0, 30000);

	dssdInfo *dsf[2], *dsb[2];
	dsf[0] = dsf[1] = dsb[0] = dsb[1] = nullptr;
	Double_t ft[2];
	for (auto ift = mapfs.begin(); ift != mapfs.end(); ++ift) {

		// front loop
		int fhit = ift->second.Hit;

		if (fhit > 2) continue;
		// fhit == 1 or 2
		dsf[0] = &(ift->second);
		ft[0] = ift->first;
		dsfb.Time = ft[0];
		dsfb.FStrip = dsf[0]->Strip;
		Double_t fse = dsf[0]->Energy;		// front sum energy
		if (fhit == 2) {					// front hit == 2
			++ift;
			dsf[1] = &(ift->second);
			ft[1] = ift->first;
			fse += dsf[1]->Energy;

			// sort by energy
			if (dsf[1]->Energy > dsf[0]->Energy) {
				dsf[1] = dsf[0];
				dsf[0] = &(ift->second);
				ft[1] = ft[0];
				ft[0] = ift->first;
				dsfb.Time = ft[0];
				dsfb.FStrip = dsf[0]->Strip;
			}
		}

		// back loop
		auto ibt = mapbs.lower_bound(ft[0]-toff-tw);
		for (; ibt != mapbs.end(); ++ibt) {
			if (ibt->first >= ft[0]-toff+tw) break;
			mapfbs.clear();

			dsb[0] = &(ibt->second);
			dsfb.BStrip = dsb[0]->Strip;
			Double_t bse = dsb[0]->Energy;

			if (dsb[0]->Hit == 1) {						// bhit == 1
				hxyse->Fill(fse-bse, fse);
				if (abs(fse - bse) < eThre) {				// energy correlation
					if (fhit == 1) {										// f1b1t1
						dsfb.FEnergy = fse;
						dsfb.BEnergy = bse;
						Double_t maxe = fse > bse ? fse : bse;
						mapfbs.insert(make_pair(maxe, dsfb));
						hxyec->Fill(fse, bse);
					} else {												// f2b1t1
						if (abs(dsf[0]->Strip - dsf[1]->Strip) == 1) {				// front pair
							dsfb.FStrip += (dsf[1]->Strip-dsf[0]->Strip) * dsf[1]->Energy / fse;
							dsfb.FEnergy = fse;
							dsfb.BEnergy = bse;
							Double_t maxe = fse > bse ? fse : bse;
							mapfbs.insert(make_pair(maxe, dsfb));
							hxyec->Fill(fse, bse);
						}
					}
				} else {
					continue;								// different energy, jump
				}
			} else {									// bhit == 2
				++ibt;
				dsb[1] = &(ibt->second);
				bse += dsb[1]->Energy;

				// sort by energy
				if (dsb[1]->Energy > dsb[0]->Energy) {
					dsb[1] = dsb[0];
					dsb[0] = &(ibt->second);
					dsfb.BStrip = dsb[0]->Strip;
				}

				hxyse->Fill(fse-bse, fse);

				if (abs(fse - bse) < eThre) {				// energy correlation
					// whether is pari strips, but 2 particles, f2b2t2
					if (fhit == 2 && abs(dsf[0]->Energy-dsb[0]->Energy) < eThre && abs(dsf[1]->Energy-dsb[1]->Energy) < eThre) {
						for (int i = 0; i != 2; ++i) {
							dsfb.Time = ft[0];
							dsfb.FStrip = dsf[i]->Strip;
							dsfb.BStrip = dsb[i]->Strip;
							dsfb.FEnergy = dsf[i]->Energy;
							dsfb.BEnergy = dsb[i]->Energy;
							Double_t maxe = dsf[i]->Energy > dsb[i]->Energy ? dsf[i]->Energy : dsb[i]->Energy;
							mapfbs.insert(make_pair(maxe, dsfb));
							hxyec->Fill(dsf[i]->Energy, dsb[i]->Energy);
						}
					} else if (abs(dsb[0]->Strip - dsb[1]->Strip) == 1) {					// back pair
						if (fhit == 1) {									// f1b2t1
							dsfb.BStrip += (dsb[1]->Strip - dsb[0]->Strip) * dsb[1]->Energy / bse;
							dsfb.FEnergy = fse;
							dsfb.BEnergy = bse;
							Double_t maxe = fse > bse ? fse : bse;
							mapfbs.insert(make_pair(maxe, dsfb));
							hxyec->Fill(fse, bse);
						}
						if (fhit == 2 && abs(dsf[0]->Strip - dsf[1]->Strip) == 1) {			// f2b2t1
							// one particle, but both front and back pair strips
							dsfb.BStrip += (dsb[1]->Strip - dsb[0]->Strip) * dsb[1]->Energy / bse;
							dsfb.FEnergy = fse;
							dsfb.BEnergy = bse;
							Double_t maxe = fse > bse ? fse : bse;
							mapfbs.insert(make_pair(maxe, dsfb));
							hxyec->Fill(fse, bse);
						}
					}
				} else {									// energy different, jump
					continue;
				}

			}

			chit = 0;
			for (auto ifbt = mapfbs.begin(); ifbt != mapfbs.end(); ++ifbt) {
				cEnergy[chit] = ifbt->first;
				cTime[chit] = ifbt->second.Time;
				cfStrip[chit] = ifbt->second.FStrip;
				cbStrip[chit] = ifbt->second.BStrip;
				cfe[chit] = ifbt->second.FEnergy;
				cbe[chit] = ifbt->second.BEnergy;
				++chit;
			}
			opt->Fill();
		}
	}




// 		if (fhit == 1) {					// front hit == 2, but maybe contains pair strips
// 			dsf[0] = &(ift->second);
// 			ft[0] = ift->first;
// 			++ift;
// 			dsf[1] = &(ift->second);
// 			ft[1] = ift->first;

// 			if (abs(dfs[0]->Strip - dsf[1]->Strip) == 1) {
// 				hit == 1;
// 				fse = dsf[0]->Energy + dsf[1]->Energy;
// 			}
// 			dsxy.Time = ft[0];
// 			dsxy.XStrip = dsf[0]->Strip + (dsf[1]->Strip-dsf[0]->Strip) * dsf[1]->Energy / fse;
// 		} else {						// front hit == 1
// 			fhit = 1;
// 			fse = ift->second.Energy;
// 			dsxy.Time = ift->first;
// 			dsxy.XStrip = ift->second.Strip;
// 			dsf[0] = &(ift->second);
// 			ft[0] = ift->first;
// 		}
// printf("%d\n", __LINE__);

// 		if (hit == 2 && abs(dsf[0]->Strip - dsf[1]->Strip) == 1) {			// confirm fhit == 1, pair strips
// 			fhit = 1;
// 			fse = dsf[0]->Energy + dsf[1]->Energy;
// 		}
// 		else fhit = 2;											// confirm fhit == 2, two particles
// printf("%d\n", __LINE__);

// 		if (fhit == 1) {										// fhit == 1
// 			auto ibt = mapbs.lower_bound(ft[0]-toff-tw);
// 			for (; ibt != mapbs.end(); ++ibt) {
// 				mapfbs.clear();
// 				if (ibt->first >= ft[0]-toff+tw) break;
// 				if (ibt->second.Hit == 1) {							// back hit == 1
// 					hxyse->Fill(fse, ibt->second.Energy);
// 					if (abs(fse-ibt->second.Energy) < eThre) {
// 						Double_t maxe = fse > ibt->second.Energy ? fse : ibt->second.Energy;
// 						dsxy.YStrip = ibt->second.Strip;
// 						dsxy.Hit = 1;
// 						dsxy.Index = 0;
// 						mapfbs.insert(make_pair(maxe, dsxy));
// 						hxyec->Fill(fse, ibt->second.Energy);
// 					}
// 				} else if (ibt->second.Hit == 2) {					// back hit == 2
// 					dsb[0] = &(ibt->second);
// 					++ibt;
// 					dsb[1] = &(ibt->second);
// 					if (abs(dsb[0]->Strip-dsb[1]->Strip) == 1) {				// pair strips
// 						Double_t bse = dsb[0]->Energy + dsb[1]->Energy;		// back sum energy
// 						hxyse->Fill(fse, bse);
// 						if (abs(fse - bse) < eThre) {
// 							Double_t maxe = fse > bse ? fse : bse;
// 							dsxy.YStrip = dsb[0]->Strip + (dsb[1]->Strip - dsb[0]->Strip) * dsb[1]->Energy / bse;
// 							dsxy.Hit = 1;
// 							dsxy.Index = 0;
// 							mapfbs.insert(make_pair(maxe, dsxy));
// 							hxyec->Fill(fse, bse);
// 						}
// 					}
// 				} else {			// back hit > 2
// 					continue;
// 				}
// 	printf("%d\n", __LINE__);

// 				chit = 0;
// 				for (auto ifbt = mapfbs.begin(); ifbt != mapfbs.end(); ++ifbt) {
// 					cEnergy[chit] = ifbt->first;
// 					cTime[chit] = ifbt->second.Time;
// 					cfStrip[chit] = ifbt->second.XStrip;
// 					cbStrip[chit] = ifbt->second.YStrip;
// 					++chit;
// 				}
// 				opt->Fill();
// 			}
// 	printf("%d\n", __LINE__);

// 		} else {											// fhit == 2, 2 particles
// 			for (int i = 0; i != 2; ++i) {
// 				auto ibt = mapbs.lower_bound(ft[i]-toff-tw);
// 				for (; ibt != mapbs.end(); ++ibt) {
// 					mapfbs.clear();
// 					if (ibt->first >= ft[i]-toff+tw) break;
// 					if (ibt->second.Hit == 1) {				// back hit == 1, no correlation
// 						continue;
// 					} else if (ibt->second.Hit == 2) {		// back hit == 2, assume 2 particle
// 						dsb[0] = &(ibt->second);
// 						++ibt;
// 						dsb[1] = &(ibt->second);
// 						hxyse->Fill(dsf[0]->Energy, dsb[0]->Energy);
// 						hxyse->Fill(dsf[1]->Energy, dsb[1]->Energy);
// 						hxyse->Fill(dsf[0]->Energy, dsb[1]->Energy);
// 						hxyse->Fill(dsf[1]->Energy, dsb[0]->Energy);
// 						if (abs(dsb[0]->Energy-dsf[0]->Energy)<eThre && abs(dsb[1]->Energy-dsf[1]->Energy)<eThre) {
// 							Double_t maxe = dsb[0]->Energy > dsf[0]->Energy ? dsb[0]->Energy : dsf[0]->Energy;
// 							dsxy.Time = ft[0];
// 							dsxy.XStrip = dsf[0]->Strip;
// 							dsxy.YStrip = dsb[0]->Strip;
// 							dsxy.Hit = 2;
// 							dsxy.Index = 0;
// 							mapfbs.insert(make_pair(maxe, dsxy));

// 							maxe = dsb[1]->Energy > dsf[1]->Energy ? dsb[1]->Energy : dsf[1]->Energy;
// 							dsxy.Time = ft[1];
// 							dsxy.XStrip = dsf[1]->Strip;
// 							dsxy.YStrip = dsb[1]->Strip;
// 							dsxy.Hit = 2;
// 							dsxy.Index = 1;
// 							mapfbs.insert(make_pair(maxe, dsxy));

// 							hxyec->Fill(dsf[0]->Energy, dsb[0]->Energy);
// 							hxyec->Fill(dsf[1]->Energy, dsb[1]->Energy);
// 						} else if (abs(dsb[0]->Energy-dsf[1]->Energy)<eThre && abs(dsb[1]->Energy-dsf[0]->Energy)<eThre) {
// 							Double_t maxe = dsb[0]->Energy > dsf[1]->Energy ? dsb[0]->Energy : dsf[1]->Energy;
// 							dsxy.Time = ft[1];
// 							dsxy.XStrip = dsf[1]->Strip;
// 							dsxy.YStrip = dsb[0]->Strip;
// 							dsxy.Hit = 2;
// 							dsxy.Index = 0;
// 							mapfbs.insert(make_pair(maxe, dsxy));

// 							maxe = dsb[1]->Energy > dsf[0]->Energy ? dsb[1]->Energy : dsf[0]->Energy;
// 							dsxy.Time = ft[0];
// 							dsxy.XStrip = dsf[0]->Strip;
// 							dsxy.YStrip = dsb[1]->Strip;
// 							dsxy.Hit = 2;
// 							dsxy.Index = 1;
// 							mapfbs.insert(make_pair(maxe, dsxy));

// 							hxyec->Fill(dsf[0]->Energy, dsb[1]->Energy);
// 							hxyec->Fill(dsf[1]->Energy, dsb[0]->Energy);
// 						}
// 					} else {					// hit > 2
// 						continue;
// 					}

// 					chit = 0;
// 					for (auto ifbt = mapfbs.begin(); ifbt != mapfbs.end(); ++ifbt) {
// 						cEnergy[chit] = ifbt->first;
// 						cTime[chit] = ifbt->second.Time;
// 						cfStrip[chit] = ifbt->second.XStrip;
// 						cbStrip[chit] = ifbt->second.YStrip;
// 						++chit;
// 					}
// 					opt->Fill();
// 				}
// 			}
// 		}
// 	}

	hxyse->Write();
	hxyec->Write();
	opt->Write();
	return 0;
}

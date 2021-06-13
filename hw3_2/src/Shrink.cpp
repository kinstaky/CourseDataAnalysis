#include "Shrink.h"

using std::multimap;
using std::pair;
using std::vector;
using std::make_pair;

DSSDShrink::DSSDShrink(TTree *it, TTree *st, TTree *rt) {
	ipt = it;
	stree = st;
	rtree = rt;
	dssdNum = 0;
}

DSSDShrink::~DSSDShrink() {

}

int DSSDShrink::AddDSSD(const char *nn, int xs, int ys, int thre) {
	names.push_back(nn);
	xStrip.push_back(xs);
	yStrip.push_back(ys);
	threshold.push_back(thre);
	dssdNum++;
	return 0;
}


int DSSDShrink::Shrink() {
	std::vector<int> sumXStrip, sumYStrip;
	sumXStrip.push_back(0);
	sumYStrip.push_back(0);
	for (int i = 1; i != dssdNum+1; ++i) {
		sumXStrip.push_back(sumXStrip[i-1] + xStrip[i-1]);
		sumYStrip.push_back(sumYStrip[i-1] + yStrip[i-1]);
	}

	// shrinking tree
	vector<Int_t> nhit;
	vector<Double_t*> nxe, nxs, nye, nys;
	// new data
	vector<Double_t*> me;			// max(xe, ye)
	vector<Int_t*> sFlag;			// xhit-yhit
									// 0x11: x-1, y-1, shrink-1
									// 0x21: x-2, y-1, shrink-1
									// 0x12: x-1, y-2, shrink-1
									// 0x22: x-2, y-2, shrink-1
	// input tree and rest tree
	vector<Int_t> xhit, yhit;
	vector<Int_t*> xs, ys;
	vector<Double_t*> xe, ye;

	for (int i = 0; i != dssdNum+1; ++i) {
		// data to fill
		nhit.push_back(0);
		nxe.push_back(new Double_t[xStrip[i]]);
		nxs.push_back(new Double_t[xStrip[i]]);

		nye.push_back(new Double_t[yStrip[i]]);
		nys.push_back(new Double_t[yStrip[i]]);

		me.push_back(new Double_t[xStrip[i]]);
		sFlag.push_back(new Int_t[xStrip[i]]);


		// data input
		xhit.push_back(0);
		xe.push_back(new Double_t[xStrip[i]]);
		xs.push_back(new Int_t[xStrip[i]]);

		yhit.push_back(0);
		ye.push_back(new Double_t[yStrip[i]]);
		ys.push_back(new Int_t[yStrip[i]]);
	}


	for (int i = 0; i != dssdNum; ++i) {
		TString hitStr = names[i]+"hit";
		TString xhitStr = names[i]+"xhit";
		TString yhitStr = names[i]+"yhit";
		TString xeStr = names[i]+"xe";
		TString yeStr = names[i]+"ye";
		TString xsStr = names[i]+"xs";
		TString ysStr = names[i]+"ys";
		TString meStr = names[i]+"me";				// max(xe, ye)
		TString flagStr = names[i]+"sFlag";			// shrinking flag

		// output shrink tree brnach
		stree->Branch(hitStr.Data(), &nhit[i], (hitStr+"/I").Data());
		stree->Branch(xeStr.Data(), nxe[i], (xeStr+"["+hitStr+"]/D").Data());
		stree->Branch(xsStr.Data(), nxs[i], (xsStr+"["+hitStr+"]/D").Data());
		stree->Branch(yeStr.Data(), nye[i], (yeStr+"["+hitStr+"]/D").Data());
		stree->Branch(ysStr.Data(), nys[i], (ysStr+"["+hitStr+"]/D").Data());
		stree->Branch(meStr.Data(), me[i], (meStr+"["+hitStr+"]/D").Data());
		stree->Branch(flagStr.Data(), sFlag[i], (flagStr+"["+hitStr+"]/I").Data());

		// output rest tree branch
		rtree->Branch(xhitStr.Data(), &xhit[i], (xhitStr+"/I").Data());
		rtree->Branch(xeStr.Data(), xe[i], (xeStr+"["+xhitStr+"]/D").Data());
		rtree->Branch(xsStr.Data(), xs[i], (xsStr+"["+xhitStr+"]/I").Data());
		rtree->Branch(yhitStr.Data(), &yhit[i], (yhitStr+"/I").Data());
		rtree->Branch(yeStr.Data(), ye[i], (yeStr+"["+yhitStr+"]/D").Data());
		rtree->Branch(ysStr.Data(), ys[i], (ysStr+"["+yhitStr+"]/I").Data());

		// input tree branch
		ipt->SetBranchAddress(xhitStr.Data(), &xhit[i]);
		ipt->SetBranchAddress(xeStr.Data(), xe[i]);
		ipt->SetBranchAddress(xsStr.Data(), xs[i]);
		ipt->SetBranchAddress(yhitStr.Data(), &yhit[i]);
		ipt->SetBranchAddress(yeStr.Data(), ye[i]);
		ipt->SetBranchAddress(ysStr.Data(), ys[i]);
	}

	// // data to fill
	// Int_t *nxhit = new Int_t[dssdNum];
	// Double_t *nxe = new Double_t[sumXStrip[dssdNum]];
	// Double_t *nxs = new Double_t[sumXStrip[dssdNum]];
	// Double_t *nxres = new Double_t[sumXStrip[dssdNum]];

	// Int_t *nyhit = new Double_t[sumYStrip[dssdNum]];
	// Double_t *nye = new Double_t[sumYStrip[dssdNum]];
	// Double_t *nys = new Double_t[sumYStrip[dssdNum]];
	// Double_t *nyres = new Double_t[sumYStrip[dssdNum]];

	// // data input
	// Int_t *xhit = new Int_t[dssdNum];
	// Double_t *xe = new Double_t[sumXStrip[dssdNum]];
	// Int_t *xs = new Int_t[sumXStrip[dssdNum]];
	// Double_t *xres = new Double_t[sumXStrip[dssdNum]];

	// Int_t *yhit = new Int_t[dssdNum];
	// Double_t *ye = new Double_t[sumYStrip[dssdNum]];
	// Int_t *ys = new Int_t[sumYStrip[dssdNum]];
	// Double_t *yres = new Double_t[sumYStrip[dssdNum]];


	// // branch
	// for (int i = 0; i != dssdNum; ++i) {
	// 	TString xhitStr = names[i]+"xhit";
	// 	TString yhitStr = names[i]+"yhit";
	// 	TString xeStr = names[i]+"xe";
	// 	TString yeStr = names[i]+"ye";
	// 	TString xsStr = names[i]+"xs";
	// 	TString ysStr = names[i]+"ys";
	// 	TString xresStr = names[i]+"xres";
	// 	TString yresStr = names[i]+"yres";

	// 	// output shrink tree brnach
	// 	stree->Branch(xhitStr.Data(), nxhit+i, (xhitStr+"/I").Data());
	// 	stree->Branch(xeStr.Data(), nxe+sumXStrip[i], (xeStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(xsStr.data(), nxs+sumXStrip[i], (xsStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(xresStr.data(), nxres+sumXStrip[i], (xresStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(yhitStr.data(), nyhit+i, (yhitStr+"/I").Data());
	// 	stree->Branch(yeStr.data(), nye+sumYStrip[i], (yeStr+"["+yhitStr+"]/D").Data());
	// 	stree->Branch(ysStr.data(), nys+sumYStrip[i], (ysStr+"["+yhitStr+"]/D").Data());
	// 	stree->Branch(yresStr.data(), nyres+sumYStrip[i], (yresStr+"["+yhitStr+"]/D").Data());

	// 	// output rest tree branch
	// 	stree->Branch(xhitStr.Data(), nxhit+i, (xhitStr+"/I").Data());
	// 	stree->Branch(xeStr.Data(), nxe+sumXStrip[i], (xeStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(xsStr.data(), nxs+sumXStrip[i], (xsStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(xresStr.data(), nxres+sumXStrip[i], (xresStr+"["+xhitStr+"]/D").Data());
	// 	stree->Branch(yhitStr.data(), nyhit+i, (yhitStr+"/I").Data());
	// 	stree->Branch(yeStr.data(), nye+sumYStrip[i], (yeStr+"["+yhitStr+"]/D").Data());
	// 	stree->Branch(ysStr.data(), nys+sumYStrip[i], (ysStr+"["+yhitStr+"]/D").Data());
	// 	stree->Branch(yresStr.data(), nyres+sumYStrip[i], (yresStr+"["+yhitStr+"]/D").Data());

	// 	// input tree branch
	// 	ipt->SetBranchAddress(xhitStr.Data(), xhit+i);
	// 	ipt->SetBranchAddress(xeStr.Data(), xe+sumXStrip[i]);
	// 	ipt->SetBranchAddress(xsStr.Data(), xs+sumXStrip[i]);
	// 	ipt->SetBranchAddress(xresStr.Data(), xres+sumXStrip[i]);
	// 	ipt->SetBranchAddress(yhitStr.Data(), yhit+i);
	// 	ipt->SetBranchAddress(yeStr.Data(), ye+sumYStrip[i]);
	// 	ipt->SetBranchAddress(ysStr.Data(), ys+sumYStrip[i]);
	// 	ipt->SetBranchAddress(yresStr.Data(), yres+sumYStrip[i]);
	// }

	multimap<Double_t, Event> nEvent;
	vector<Double_t> txe, tye;
	vector<Int_t> txs, tys;

	Long64_t nentry = ipt->GetEntries();
	Long64_t nentry100 = nentry / 100;
	printf("shrinking   0%%");
	fflush(stdout);
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		ipt->GetEntry(jentry);
		bool valid = true;

		// printf("===========\nentry %lld\n===========\n", jentry);
		for (int i = 0; i != dssdNum; ++i) {
			nEvent.clear();
			Int_t nxhit = xhit[i];
			Int_t nyhit = yhit[i];

			// temporay data to be convenient to opearte
			txe.clear();
			tye.clear();
			txs.clear();
			tys.clear();
			for (int j = 0; j != nxhit; ++j) {
				txe.push_back(xe[i][j]);
				txs.push_back(xs[i][j]);
			}
			for (int j = 0; j != nyhit; ++j) {
				tye.push_back(ye[i][j]);
				tys.push_back(ys[i][j]);
			}

			while (nxhit > 0 && nyhit > 0 && valid) {
				// printf("loop nxhit %d, nyhit %d\n", nxhit, nyhit);
				int findX = -1;
				int findY = -1;
				int findY2 = -1;						// 找到相邻条
				Double_t ne = txe[0];
				for (int j = 1; j < nxhit; ++j) {
					// printf("finding x, j = %d\n", j);
					if (abs(txs[0]-txs[j]) == 1) {
						if (findX >= 0) {
							// printf("findX %d, next 3, invalid, break\n", xs[i][j]);
							valid = false;
							break;
						} else {
							findX = j;
							ne += txe[j];
							// printf("find first x %d\n", findX);
						}
					}
				}

				if (!valid) break;

				Double_t diffe = 1000.0;
				// 寻找y面单条对应
				for (int j = 0; j < nyhit; ++j) {
					// printf("finding y j = %d\n", j);
					Double_t de = fabs(tye[j] - ne);
					if (de < threshold[i]) {
						if (de < diffe) {					// 无论是不是找到的第一条都可以以此判断为能量差最小的合适条
							// printf("find single y[%d], de %lf\n", ys[i][j], de);
							findY = j;
							diffe = de;
						}
					} else if (tye[j] < ne) {
						break;								// 能量已经小于ne，不用再往下找了
					}
				}
				// 寻找y面对应的相邻条
				for (int j = 0; j < nyhit-1; ++j) {
					// printf("finding y1 j = %d\n", j);
					for (int k = j+1; k != nyhit; ++k) {
						// printf("finding y2 k = %d\n", k);
						Double_t de = fabs(tye[j]+tye[k]-ne);
						if (de < threshold[i] && de < diffe) {
							// printf("found pair strips, y[%d] and y[%d], de %lf\n", ys[i][j], ys[i][k], de);
							findY = j;
							findY2 = k;
							diffe = de;
						}
					}
					if (tye[j]*2.0+120.0 < ne) {
						break;							// 能量已经小于ne，不用再往下找了
					}
				}

				if (findX > 0 && findY2 > 0 && fabs(txe[0]-tye[findY]) < 10 && fabs(txe[findX]-tye[findY2]) < 10) {
					// printf("find 222\n");
					Event event;
					event.Xe = txe[findX];
					event.Ye = tye[findY2];
					event.XStrip = txs[findX];
					event.YStrip = tys[findY2];
					event.Flag = 0x11;
					Double_t maxe = event.Xe > event.Ye ? event.Xe : event.Ye;
					nEvent.insert(make_pair(maxe, event));
					// printf("insert %lf, xe %lf, xs %lf, ye %lf, ys %lf, flag %x\n", maxe, event.Xe, event.XStrip, event.Ye, event.YStrip, event.Flag);

					// moving data, erase choosed data
					txe.erase(txe.begin()+findX);
					txs.erase(txs.begin()+findX);
					tye.erase(tye.begin()+findY2);
					tys.erase(tys.begin()+findY2);

					// turn back to single strip，化归到只有单条对应
					nxhit -= 1;
					nyhit -= 1;
					findX = -1;
					findY2 = -1;
				}

				Event event;
				// printf("xe[i][0] = %lf, xe[i][findX] = %lf\n", xe[i][0], findX == -1 ? 0.0 : xe[i][findX]);
				// 单条
				if (findX == -1) {
					event.Xe = txe[0];
					event.XStrip = txs[0];
					event.Flag |= 0x10;
					// erase x data
					txe.erase(txe.begin());
					txs.erase(txs.begin());

					nxhit -= 1;
				} else {
					event.Xe = txe[0] + txe[findX];
					event.XStrip = txs[0] + txe[findX] / event.Xe * (txs[findX] - txs[0]);
					event.Flag |= 0x20;
					// erase x data
					txe.erase(txe.begin()+findX);
					txs.erase(txs.begin()+findX);
					txe.erase(txe.begin());
					txs.erase(txs.begin());

					nxhit -= 2;
				}

				if (findY == -1 && findY2 == -1) {		// 没有找到对应的能量
					valid = false;
				} else if (findY2 == -1) {				// 单条对应
					event.Ye = tye[findY];
					event.YStrip = tys[findY];
					event.Flag |= 0x01;
					// erase y data
					tye.erase(tye.begin()+findY);
					tys.erase(tys.begin()+findY);

					nyhit -= 1;
				} else {
					event.Ye = tye[findY] + tye[findY2];
					event.YStrip = tys[findY] + tye[findY2] / event.Ye * (tys[findY2] - tys[findY]);
					event.Flag |= 0x02;
					// erase y data
					tye.erase(tye.begin()+findY2);
					tys.erase(tys.begin()+findY2);
					tye.erase(tye.begin()+findY);
					tys.erase(tys.begin()+findY);

					nyhit -= 2;
				}
				Double_t maxe = event.Xe > event.Ye ? event.Xe : event.Ye;
				nEvent.insert(make_pair(maxe, event));
				// printf("insert %lf, xe %lf, xs %lf, ye %lf, ys %lf, flag %x\n", maxe, event.Xe, event.XStrip, event.Ye, event.YStrip, event.Flag);
			}


			if (nxhit != nyhit) {
				// printf("invalid nxhit = %d, nyhit = %d\n", nxhit, nyhit);
				if (nxhit == 0) {
					for (int j = 0; j != nyhit; ++j) {
						if (tye[j] < 100 + threshold[i]) {
							Event event;
							event.Xe = 0.0;
							event.Ye = tye[j];
							event.XStrip = -1;
							event.YStrip = tys[j];
							event.Flag = 0x01;
							nEvent.insert(make_pair(event.Ye, event));
						} else {
							valid = false;
							break;
						}
					}
				} else if (nyhit == 0) {
					for (int j = 0; j != nxhit; ++j) {
						if (txe[j] < 100 + threshold[i]) {
							Event event;
							event.Xe = txe[j];
							event.Ye = 0.0;
							event.XStrip = txs[j];
							event.YStrip = -1;
							event.Flag = 0x10;
							nEvent.insert(make_pair(event.Xe, event));
						} else {
							valid = false;
							break;
						}
					}
				} else {
					valid = false;
				}
			}

			if (valid) {
				nhit[i] = nEvent.size();
				int index = 0;
				for (auto it = nEvent.rbegin(); it != nEvent.rend(); ++it) {
					Event &ev = it->second;
					me[i][index] = it->first;
					nxe[i][index] = ev.Xe;
					nye[i][index] = ev.Ye;
					nxs[i][index] = ev.XStrip;
					nys[i][index] = ev.YStrip;
					sFlag[i][index] = ev.Flag;
					++index;
					// printf("fill me %lf, xe %lf, xs %lf, ye %lf, ys %lf, flag %x\n", it->first, ev.Xe, ev.XStrip, ev.Ye, ev.YStrip, ev.Flag);
				}
			}
		}

		if (valid) {
			stree->Fill();
		} else {
			rtree->Fill();
		}

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", (jentry+1)/nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");


	for (int i = 0; i != dssdNum; ++i) {
		delete[] nxe[i];
		delete[] nxs[i];
		delete[] nye[i];
		delete[] nys[i];
		delete[] xs[i];
		delete[] ys[i];
		delete[] xe[i];
		delete[] ye[i];
		delete[] me[i];
		delete[] sFlag[i];
	}
	return 0;
}

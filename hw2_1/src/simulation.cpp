#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TString.h>
#include <ctime>

// detector data
const int DETECTORS = 5;
const int PARTICLES = 1000000;
const Double_t CHI_THRESHOLD = 10.0 * (Double_t(DETECTORS) - 2.0);
// position
const Double_t xz[DETECTORS]{-1750.7, -1713.3, -1250.7, -1213.3, -381.6};
const Double_t yz[DETECTORS]{-1742.1, -1721.9, -1242.1, -1221.9, -373.4};
// resolution
const Double_t sx[DETECTORS]{0.336971, 0.365249, 0.399483, 0.405394, 0.202764};
const Double_t sy[DETECTORS]{0.365336, 0.375106, 0.42539, 0.415491, 0.234976};
// step
const int MAXRUN = 100;
const double eps = 1e-6;
const int MAXBATCH = 10;

// data
Double_t rx[PARTICLES][DETECTORS];	// “真实”位置
// global variables
TRandom3 *gr;

void genTrace(TH1D *hNear, TH1D *hFar, Double_t rp[][DETECTORS], const Double_t *z) {
	for (int i = 0; i != PARTICLES; ++i) {
		Double_t pFar = hFar->GetRandom(gr);
		Double_t pNear = hNear->GetRandom(gr);
		Double_t k = (pNear - pFar ) / 1800.0;
		Double_t b = pNear;
		for (int j = 0; j != DETECTORS; ++j) {
			rp[i][j] = k * z[j] + b;			// “真实”位置
		}
	}
}


Double_t sFit(Double_t *fx, Double_t *fy, Double_t &k, Double_t &b) {
	const Int_t n = DETECTORS;
	Double_t sumx = 0.0;
	Double_t sumy = 0.0;
	Double_t sumxy = 0.0;
	Double_t sumx2 = 0.0;
	for (Int_t i = 0; i != n; ++i) {
		Double_t x = fx[i];
		Double_t y = fy[i];
		sumx += x;
		sumy += y;
		sumxy += x * y;
		sumx2 += x * x;
	}
	Double_t dn = Double_t(n);
	k = (sumxy - sumx*sumy/dn) / (sumx2 - sumx*sumx/dn);
	b = (sumy - k*sumx) / dn;
	Double_t chi2 = 0.0;
	for (Int_t i = 0; i != n; ++i) {
		Double_t x = fx[i];
		Double_t y = fy[i];
		Double_t t = y - k*x - b;
		chi2 += t * t;
	}
	return chi2;
}


// 模拟粒子打在探测器上，因为分辨位置偏移，利用偏移后的位置重建径迹，返回残差分布
void hitAndTrace(TH1D **hdp, Double_t rp[][DETECTORS], Double_t *ss, const Double_t *z) {
	for (int j = 0; j != DETECTORS; ++j) {
		TString name = Form("hdp%d", j);
		hdp[j] = new TH1D(name.Data(), name.Data(), 200, -1, 1);
	}

	// hit and trace
	for (int i = 0; i != PARTICLES; ++i) {
		Double_t vp[DETECTORS];
		for (int j = 0; j != DETECTORS; ++j) {
			vp[j] = gr->Gaus(rp[i][j], ss[j]);			// hit and bring resolution
		}
		Double_t k, b;
		sFit((Double_t*)z, (Double_t*)vp, k, b);
		// if (sFit((Double_t*)z, (Double_t*)vp, k, b) >= CHI_THRESHOLD) continue;			// no need, the chi2 is very small
		for (int j = 0; j != DETECTORS; ++j) {
			Double_t tp = b + k * z[j];			// 用径迹拟合的位置
			hdp[j]->Fill(vp[j]-tp);					// 填充残差
		}
	}
}

void fitResidual(TH1D **hdp, Double_t *ns) {
	TF1 *f1 = new TF1("f1", "gaus", -0.6, 0.6);
	for (int j = 0; j != DETECTORS; ++j) {
		hdp[j]->Fit(f1, "QRN");
		ns[j] = f1->GetParameter(2);
	}
	delete f1;
}

typedef struct stepParameter {
	Double_t alpha;				// learning rate ??
	Double_t loss;				// loss ??
} stepParameter;
stepParameter sp;
// 根据拟合得到的新的残差调整分辨率，ns - 当前残差sigma， ts - 目标残差sigma， sr - 模拟的分辨率
void step(Double_t *ns, const Double_t *ts, Double_t *sr) {
	for (int j = 0; j != DETECTORS; ++j) {
		ns[j] = ts[j] - ns[j];
		sr[j] += sp.alpha * ns[j];
	}
	return;
}
// 计算loss函数，ns - 当前残差sigma，ts - 目标残差sigma
Double_t lossFunc(Double_t *ns, const Double_t *ts) {
	Double_t l = 0.0;
	for (int j = 0; j != DETECTORS; ++j) {
		Double_t t = ns[j] - ts[j];
		l += t * t;
	}
	if (l > sp.loss) {
		sp.alpha *= 0.8;
	}
	sp.loss = l;
	return l;
}


int simulation() {
	// output
	Double_t minLoss = 1.0;					// minimum loss
	Double_t minSsx[DETECTORS];				// ssx for minimum loss
	TFile *opf = new TFile("../data/simRecord.root", "recreate");
	TTree *opt = new TTree("tree", "tree for record simulation");
	opt->Branch("loss", &minLoss, "loss/D");
	opt->Branch("res", &minSsx, TString::Format("res[%d]/D", DETECTORS).Data());


	// random engine
	gr = new TRandom3(0);


	int batchN = 0;
	while (batchN < MAXBATCH) {
		TFile *ipf = new TFile("../data/tracking.root");
		TTree *ipt = (TTree*)ipf->Get("tree");


		// 获取分布
		ipt->Draw("bx>>hxn(120, -30, 30", "c2nx < 10", "goff");
		TH1D *hxn = (TH1D*)gDirectory->Get("hxn");
		ipt->Draw("-1800*kx+bx>>hxf(120, -40, 40)", "c2nx < 10", "goff");
		TH1D *hxf = (TH1D*)gDirectory->Get("hxf");


		clock_t t = clock();
		// 生成径迹和计算径迹和探测器的交点
		genTrace(hxn, hxf, rx, xz);

		ipf->Close();


		Double_t ssx[DETECTORS];				// present simulating resolution
		// semi-random init ssx
		for (int j = 0; j != DETECTORS; ++j) {
			ssx[j] = gr->Gaus(sx[j], 0.1);
		}
		TH1D *hdx[DETECTORS];				// residual histogram
		TH1D *mhdx[DETECTORS];				// residual histogram for minium loss
		Double_t nsx[DETECTORS];			// new residual sigma

		sp.alpha = 1.0;
		sp.loss = minLoss;
		int run = 0;
		while (true) {
			hitAndTrace((TH1D**)hdx, rx, ssx, xz);		// simulate particles hit the detetors and trace it
			fitResidual((TH1D**)hdx, nsx);				// fit and get resolution
			Double_t loss = lossFunc(nsx, sx);			// calculate loss
			if (run == 0) {
				for (int j = 0; j != DETECTORS; ++j) {
					mhdx[j] = hdx[j];
					mhdx[j]->SetName(TString::Format("hdx%d_%03d", j, batchN).Data());
					minSsx[j] = ssx[j];
				}
				minLoss = loss;
			} else if (loss < minLoss) {
				minLoss = loss;
				for (int j = 0; j != DETECTORS; ++j) {
					delete mhdx[j];
					mhdx[j] = hdx[j];
					mhdx[j]->SetName(TString::Format("hdx%d_%03d", j, batchN).Data());
					minSsx[j] = ssx[j];
				}
			} else {
				for (int j = 0; j != DETECTORS; ++j) {		// delete hdx
					delete hdx[j];
				}
			}
			if (run == MAXRUN || loss < eps) break;
			step(nsx, sx, ssx);							// automatically change the simulating resolution
			printf("batch %d, run %d, loss %lf, rate %lf\n", batchN, run, loss, sp.alpha);		// show the process
			run++;
		}


		printf("%f s\n", float(clock()-t) / CLOCKS_PER_SEC);

		printf("\ntarget    present   resolution\n");
		for (int j = 0; j != DETECTORS; ++j) {
			printf("%lf  %lf  %lf\n", sx[j], nsx[j], ssx[j]);
		}

		opf->cd();
		for (int j = 0; j != DETECTORS; ++j) {
			mhdx[j]->Write();
			delete mhdx[j];
		}
		opt->Fill();


		++batchN;
	}

	opt->Write();
	opf->Close();

	// tree->Draw("by>>hyn(120, -30, 30)", "c2ny < 10", "goff");
	// tree->Draw("-1800*ky+by>>hyf(120, -60, 60)", "c2ny < 10", "goff");


	return 0;
}

void test(Double_t *ssx) {
	gr = new TRandom3(0);

	TFile *ipf = new TFile("../data/tracking.root");
	TTree *ipt = (TTree*)ipf->Get("tree");


	// 获取分布
	ipt->Draw("bx>>hxn(120, -30, 30", "c2nx < 10", "goff");
	TH1D *hxn = (TH1D*)gDirectory->Get("hxn");
	ipt->Draw("-1800*kx+bx>>hxf(120, -40, 40)", "c2nx < 10", "goff");
	TH1D *hxf = (TH1D*)gDirectory->Get("hxf");


	clock_t t = clock();
	// 生成径迹和计算径迹和探测器的交点
	genTrace(hxn, hxf, rx, xz);

	ipf->Close();

	TH1D *hdx[DETECTORS];				// residual histogram
	Double_t nsx[DETECTORS];			// new residual sigma

	hitAndTrace((TH1D**)hdx, rx, ssx, xz);		// simulate particles hit the detetors and trace it
	fitResidual((TH1D**)hdx, nsx);				// fit and get resolution
	Double_t loss = lossFunc(nsx, sx);			// calculate loss
	for (int j = 0; j != DETECTORS; ++j) {		// delete hdx
		delete hdx[j];
	}

	// show the result
	printf("test %d, loss %lf\n", 0, loss);
	printf("%f s\n", float(clock()-t) / CLOCKS_PER_SEC);
	printf("\n| target   | present  | resolution |\n");
	printf("| -------- | -------- | -----------|\n");
	for (int j = 0; j != DETECTORS; ++j) {
		printf("| %lf | %lf | %lf   |\n", sx[j], nsx[j], ssx[j]);
	}
}



#ifndef __CINT__
int main(int argc, char **argv) {
	// simulation();
	Double_t tsx[5]{0.388706, 0.472974, 0.420475, 0.421027, 0.762936};
	test(tsx);
	return 0;
}
#endif
#include "Tracking.h"

Tracking::Tracking(TTree *tree)
	: trackingBase(tree) {
}


Tracking::~Tracking() {
}



void Tracking::setBranch(TTree *tree) {
	tree->Branch("xx", &xx, "xx[5]/D");
	tree->Branch("xz", &xz, "xz[5]/D");
	tree->Branch("yy", &yy, "yy[5]/D");
	tree->Branch("yz", &yz, "yz[5]/D");
	tree->Branch("anode", &anode, "anode[5]/D");

	tree->Branch("dx", &dx, "dx[5]/D");
	tree->Branch("dy", &dy, "dy[5]/D");

	tree->Branch("tx", &tx, "tx/D");
	tree->Branch("ty", &ty, "ty/D");
	tree->Branch("tz", &tz, "tz/D");

	tree->Branch("ptx", &ptx, "ptx/D");
	tree->Branch("pty", &pty, "pty/D");

	tree->Branch("c2nx", &c2nx, "c2nx/D");
	tree->Branch("c2ny", &c2ny, "c2ny/D");

	tree->Branch("beamTrig", &beamTrig, "beamTrig/I");
	tree->Branch("must2Trig", &must2Trig, "must2Trig/I");

	tree->Branch("targetX", &targetX, "targetX");
	tree->Branch("targetY", &targetY, "targetY");

	tree->Branch("numX", &numX, "numX/I");
	tree->Branch("numY", &numY, "numY/I");

	tree->Branch("trackX", &trackX, "trackX/I");
	tree->Branch("trackY", &trackY, "trackY/I");

	tree->Branch("bx", &bx, "bx/D");
	tree->Branch("kx", &kx, "kx/D");
	tree->Branch("by", &by, "by/D");
	tree->Branch("ky", &ky, "ky/D");

	tree->Branch("fxx", &fxx, "fxx[5]/D");
	tree->Branch("fyy", &fyy, "fyy[5]/D");
}


Int_t Tracking::trackInit() {
	tx = -999;
	ty = -999;

	for (int i = 0; i != 5; ++i) {
		xx[i] = PPACF8[i][0];
		xz[i] = PPACF8[i][2];
		yy[i] = PPACF8[i][1];
		yz[i] = PPACF8[i][3];
		anode[i] = PPACF8[i][4];
	}

	// 判断是否满足重建流程：
	// 	1. 遍历所有探测器，记录有信号的探测器总数量和组合
	//	2. 如果数量等于1，不要
	//	3. 如果数量为2，且所有信号的都在PPAC1或者PPAC2里面，不要
	//	4. 其他都要

	// check x
	numX = 0;
	trackX = 0;
	for (int i = 0; i != 4; ++i) {
		if (abs(xx[i]) < 120) {
			numX += 1;
			trackX |= 1 << i;
		}
	}
	// PPAC3
	if (abs(xx[4]) < 50) {
		numX += 1;
		trackX |= 1 << 4;
	}

	if (numX == 1) numX = -1;
	if (numX == 2) {
		if (trackX & 0x3) numX = -2;			// 信号都在PPAC1
		if (trackX & 0xC) numX = -2;			// 信号都在PPAC2
	}


	// check y
	numY = 0;
	trackY = 0;
	for (int i = 0; i != 4; ++i) {
		if (abs(yy[i]) < 75) {
			numY += 1;
			trackY |= 1 << i;
		}
	}
	// PPAC3
	if (abs(yy[4]) < 50) {
		numY += 1;
		trackY |= 1 << 4;
	}

	if (numY == 1) numY = -1;
	if (numY == 2) {
		if (trackY & 0x3) numY = -2;
		if (trackY & 0xC) numY = -2;
	}


	Int_t flag = 0;
	if (numX > 0) flag |= 1;
	if (numY > 0) flag |= 2;
	return flag;
}


void Tracking::addTrace(TH2D *h, Double_t k, Double_t b, Int_t minBin, Int_t maxBin) {
	if (h == nullptr) return;
	if (minBin >= maxBin) return;

	for (int i = minBin; i != maxBin; i+=10) {
		h->Fill(i, (Int_t)(i*k*b));
	}
	return;
}

void Tracking::Loop() {
	TH2D *htf8xz = new TH2D("htf8xz", "xz track by ppac", 220, -2000, 200, 300, -150, 150);
	TH2D *htf8yz = new TH2D("htf8yz", "yz track by ppac", 220, -2000, 200, 300, -150, 150);

	TFile *opf = new TFile("../data/tracking.root", "recreate");
	TTree *tree = new TTree("tree", "ppac tracking");
	setBranch(tree);


	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry = 0; jentry != nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		// x and y tracked
		if (trackInit() != 3) continue;

		TGraph *grx = new TGraph;
		int points = 0;
		for (int i = 0; i != 5; ++i) {
			if (trackX & 1<<i) {
				grx->SetPoint(points, xz[i], xx[i]);
				points++;
			}
		}
		TF1 *fx = new TF1("fx", "pol1", -2000, 100);
		TFitResultPtr r = grx->Fit(fx, "SQ");
		bx = fx->GetParameter(0);
		kx = fx->GetParameter(1);
		// residual and fitted position
		for (int i = 0; i != 5; ++i) {
			fxx[i] = fx->Eval(xz[i]);
			if (trackX & 1<<i) {
				dx[i] = xx[i] - fxx[i];
			}
			else {
				dx[i] = -999;
			}
		}
		c2nx = r->Chi2() / r->Ndf();
		delete grx;
		delete fx;

		TGraph* gry = new TGraph;
		points = 0;
		for (int i = 0; i != 5; ++i) {
			if (trackY & 1<<i) {
				gry->SetPoint(points, yz[i], yy[i]);
				points++;
			}
		}
		TF1 *fy = new TF1("fy", "pol1", -2000, 100);
		r = gry->Fit(fy, "SQ");
		by = fy->GetParameter(0);
		ky = fy->GetParameter(1);
		// residual
		for (int i = 0; i != 5; ++i) {
			fyy[i] = fy->Eval(yz[i]);
			if (trackY & 1<<i) {
				dy[i] = yy[i] - fyy[i];
			}
			else {
				dy[i] = -999;
			}
		}
		c2ny = r->Chi2() / r->Ndf();
		delete gry;
		delete fy;


		// target position
		tz = -bx / (kx + 1.0);
		tx = -tz;
		tx = bx;
		ty = ky * tz + by;
		// projection
		ptx = tx * 1.4142;			// tx / sqrt(2.0)
		pty = ty;

		// set trace
		addTrace(htf8xz, kx, bx, -1800, 100);
		addTrace(htf8yz, ky, by, -1800, 100);

		tree->Fill();
		if (jentry % 10000 == 0) std::cout << jentry << "/" << nentries << std::endl;
	}

	htf8xz->Write();
	htf8yz->Write();
	tree->Write();
	opf->Close();
}

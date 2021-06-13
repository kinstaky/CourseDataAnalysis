#include "DSSD.h"

#include <TString.h>
#include <TH2F.h>
#include <TCut.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFitResult.h>

using std::string;
using std::vector;

DSSD::DSSD(TTree *iipt, const char *nn, int xs, int ys) {
	ipt = iipt;
	Name = string(nn);
	XStrip = xs;
	YStrip = ys;
	LogFile = nullptr;
	GraphFile = nullptr;

	// init normalize parameter
	for (int i = 0; i != xs; ++i) {
		XNormPar.push_back(0.0);
		XNormPar.push_back(1.0);
		XNormPar.push_back(0.0);
		XNormPar.push_back(0.0);

		xNormFlag.push_back(false);

		TString resName;
		resName.Form("resX[%d]", i);
		hresX.push_back(new TH2F(resName.Data(), resName.Data(), 100, -50, 50, 1000, 0, 8000));

		gFitX.push_back(new TGraph);
	}
	for (int i = 0; i != ys; ++i) {
		YNormPar.push_back(0.0);
		YNormPar.push_back(1.0);
		YNormPar.push_back(0.0);
		YNormPar.push_back(0.0);

		yNormFlag.push_back(false);

		TString resName;
		resName.Form("resY[%d]", i);
		hresY.push_back(new TH2F(resName.Data(), resName.Data(), 100, -50, 50, 1000, 0, 8000));

		gFitY.push_back(new TGraph);
	}

	// hit1Cut = TCut((Name + "xhit == 1 && " + Name + "yhit == 1").data());
	// xhit2Cut = TCut((Name + "xhit == 2 && abs(" + Name + "xs[0]-" + Name + "xs[1]) > 1").data());
	// yhit2Cut = TCut((Name + "yhit == 2 && abs(" + Name + "ys[0]-" + Name + "ys[1]) > 1").data());
	// hit2Cut = xhit2Cut && yhit2Cut;
	// hit1or2Cut = hit1Cut || hit2Cut;
}


DSSD::~DSSD() {
}

int DSSD::AddNode(int xmin, int xmax, int ymin, int ymax, int dir) {
	// 判断所给的条是否溢出
	if (xmin >= xmax || ymin >= ymax) return -1;
	if (xmin < 0 || xmin >= XStrip) return -1;
	if (xmax < 0 || xmax > XStrip) return -1;
	if (ymin < 0 || ymin >= YStrip) return -1;
	if (ymax < 0 || ymax > YStrip) return -1;
	if (dir < 0 || dir > 1) return -1;

	if (nodes.size() == 0) {
		if (dir == DIR_X2Y) {
			if (xmax != xmin+1) return -1;
			xNormFlag[xmin] = true;
		} else {
			if (ymax != ymin+1) return -1;
			yNormFlag[ymin] = true;
		}
	}
	// 保存节点信息
	NodeInfo nf;
	nf.direction = dir;
	nf.xmin = xmin;
	nf.xmax = xmax;
	nf.ymin = ymin;
	nf.ymax = ymax;
	nodes.push_back(nf);
	return 0;
}

int DSSD::genGraph() {
	GraphFile->cd();
	TDirectory *subDir = GraphFile->GetDirectory(Name.data());
	if (!subDir) {
		subDir = new TDirectoryFile(Name.data(), Name.data());
	}
	subDir->cd();
	for (int x = 0; x != XStrip; ++x) {
		for (int y = 0; y != YStrip; ++y) {
			TString gName;
			gName.Form("gX[%d]Y[%d]", x, y);
			TGraph *g = (TGraph*)subDir->Get(gName.Data());
			if (g) continue;

			// TString varStr;
			// varStr.Form("y[%d]:x[%d]", y, x);
			// TCut xCut, yCut;
			// if (x == 0) {
			// 	xCut = "x[0] > 200 && x[0] < 8192 && x[1] < 100";
			// } else if (x < XStrip-1) {
			// 	xCut = TString::Format("x[%d] < 100 && x[%d] < 100 && x[%d] > 200 && x[%d] < 8192", x-1, x+1, x, x).Data();
			// } else {
			// 	xCut = TString::Format("x[%d] < 100 && x[%d] > 200 && x[%d] < 8192", x-1, x, x).Data();
			// }
			// if (y == 0) {
			// 	yCut = "y[1] < 100 && y[0] > 200 && y[0] < 8192";
			// } else if (y < YStrip-1) {
			// 	yCut = TString::Format("y[%d] < 100 && y[%d] < 100 && y[%d] > 200 && y[%d] < 8192", y-1, y+1, y, y).Data();
			// } else {
			// 	yCut = TString::Format("y[%d] < 100 && y[%d] > 200 && y[%d] < 8192", y-1, y, y).Data();
			// }
			// ipt->Draw(varStr.Data(), xCut && yCut, "goff");

			// 1 hit cut
			TCut hit1Cut = "xhit == 1 && yhit == 1";
			// 2 hit cut
			TCut hit2Cut = "xhit == 2 && abs(xs[0]-xs[1]) > 1 && yhit == 2 && abs(ys[0]-ys[1]) > 1";
			// 选择x面需要的条
			TCut xStripCut = TString::Format("xs==%d", x).Data();
			// 选择y面需要的条
			TCut yStripCut = TString::Format("ys==%d", y).Data();
			// 能量范围
			TCut eCut = "xe < 8192 && ye < 8192";

			ipt->Draw("ye:xe", hit2Cut && xStripCut && yStripCut && eCut, "goff");
			// add points to the graph
			Long64_t rows = ipt->GetSelectedRows();
			Double_t *pointY = ipt->GetV1();
			Double_t *pointX = ipt->GetV2();

			g = new TGraph(rows, pointX, pointY);
			g->SetMarkerStyle(7);
			g->SetLineWidth(0);

			// ipt->Draw("ye:xe", hit1Cut && xStripCut && yStripCut && eCut, "goff");
			// rows = ipt->GetSelectedRows();
			// pointY = ipt->GetV1();
			// pointX = ipt->GetV2();

			// for (Long64_t ip = 0; ip != rows; ++ip) {
			// 	g->SetPoint(g->GetN(), pointX[ip], pointY[ip]);
			// }
			g->Write(gName.Data());

			printf("Graph X=%d Y=%d\n", x, y);
		}
	}
	return 0;
}


int DSSD::NormAll() {
	// set alias to adapt to different dssd
	ipt->SetAlias("x", (Name + "x").data());
	ipt->SetAlias("xhit", (Name + "xhit").data());
	ipt->SetAlias("xe", (Name + "xe").data());
	ipt->SetAlias("xs", (Name + "xs").data());

	ipt->SetAlias("y", (Name + "y").data());
	ipt->SetAlias("yhit", (Name + "yhit").data());
	ipt->SetAlias("ye", (Name + "ye").data());
	ipt->SetAlias("ys", (Name + "ys").data());

	genGraph();

	// 逐个节点进行操作
	for (long unsigned int inode = 0; inode != nodes.size(); ++inode) {
		NodeInfo &node = nodes[inode];
		printf("\n==========\nnode %ld\n==========\n", inode);
		printf("%3s %3s %7s %7s %7s\n", "dir", "i", "p0", "p1", "p2");

		// 区分刻度x面还是y面
		if (node.direction == DIR_X2Y) {
			for (int i = node.ymin; i != node.ymax; ++i) {
				// if (yNormFlag[i]) continue;
				if (NormYStrip(node.xmin, node.xmax, i) == 0) yNormFlag[i] = true;
			}
		} else {
			for (int i = node.xmin; i != node.xmax; ++i) {
				// if (xNormFlag[i]) continue;
				if (NormXStrip(node.ymin, node.ymax, i) == 0) xNormFlag[i] = true;
			}
		}
	}

	if (LogFile) {
		LogFile->cd();

		LogFile->GetDirectory(Name.data())->cd();
		gDirectory->GetDirectory("Fit")->cd();
		for (int i = 0; i != XStrip; ++i) {
			gFitX[i]->Write(TString::Format("gX[%d]", i).Data());
		}
		for (int i = 0; i != YStrip; ++i) {
			gFitY[i]->Write(TString::Format("gY[%d]", i).Data());
		}

		LogFile->GetDirectory(Name.data())->cd();
		gDirectory->GetDirectory("Residual")->cd();
		for (auto iter : hresX) {
			iter->Write();
		}
		for (auto iter : hresY) {
			iter->Write();
		}

		LogFile->GetDirectory(Name.data())->cd();
		for (int j = 0; j != 3; ++j) {
			gxNormP[j] = new TGraph;
			gyNormP[j] = new TGraph;
		}
		for (int i = 0; i != XStrip; ++i) {
			for (int j = 0; j != 3; ++j) {
				gxNormP[j]->SetPoint(i, i, XNormPar[i*4+j]);
			}
		}
		for (int i = 0; i != XStrip; ++i) {
			for (int j = 0; j != 3; ++j) {
				gyNormP[j]->SetPoint(i, i, YNormPar[i*4+j]);
			}
		}
		for (int i = 0; i != 3; ++i) {
			gxNormP[i]->Write(TString::Format("gxNormP%d", i));
			gyNormP[i]->Write(TString::Format("gyNormP%d", i));
		}
	}

	return 0;
}


int DSSD::NormXStrip(int ymin, int ymax, int x) {
	TGraph *g = gFitX[x];
	Int_t pointNum = 0;
	for (int i = ymin; i != ymax; ++i) {
		if (!yNormFlag[i]) continue;
		// // 选择x面需要的条
		// TCut xStripCut = TString::Format("xs==%d", x).Data();
		// // 选择y面需要的条
		// TCut yStripCut = TString::Format("ys==%d", i).Data();
		// // 需要画出的表达式，给参考条乘上刻度系数
		// TString varStr, hName;
		// hName.Form("h2X[%d]Y[%d]", x, i);
		// varStr.Form("%lf+%lf*ye+%lf*ye*ye:xe>>%s(1000, 0, 8000, 1000, 0, 8000)", YNormPar[i*4], YNormPar[i*4+1], YNormPar[i*4+2], hName.Data());
		// ipt->Draw(varStr.Data(), hit1or2Cut && xStripCut && yStripCut, "goff");

		// // add points to the graph
		// Long64_t rows = ipt->GetSelectedRows();
		// Double_t *pointY = ipt->GetV1();
		// Double_t *pointX = ipt->GetV2();
		// for (Long64_t j = 0; j != rows; ++j) {
		// 	g->SetPoint(pointNum, pointX[j], pointY[j]);
		// 	++pointNum;
		// }

		// if (LogFile) {
		// 	TH2F *h2xy = (TH2F*)gDirectory->Get(hName.Data());
		// 	LogFile->cd();
		// 	LogFile->GetDirectory("TH2")->cd();
		// 	h2xy->Write();
		// }

		TString gName;
		gName.Form("gX[%d]Y[%d]", x, i);
		if (!GraphFile) {
			printf("Error no GraphFile\n");
			return -1;
		}
		TGraph *gs = (TGraph*)GraphFile->Get((Name+"/"+gName).Data());
		Int_t n = gs->GetN();
		for (Int_t pi = 0; pi != n; ++pi) {
			Double_t px = gs->GetPointX(pi);
			Double_t py = gs->GetPointY(pi);
			py = YNormPar[i*4] + YNormPar[i*4+1]*py + YNormPar[i*4+2]*py*py;
			g->SetPoint(pointNum++, px, py);
		}
	}

	if (pointNum == 0) return -1;
	// fit
	// TF1 *f1 = new TF1("f1", "pol2", 200, 8000);
	// f1->SetParameter(0, 0.0);
	// f1->SetParameter(1, 1.0);
	// f1->SetParameter(2, 0.0);
	// TFitResultPtr fr = g->Fit(f1, "SQ rob=0.7+");
	// g->Fit(f1, "rob=0.7+");
	// g->SetMarkerStyle(7);
	// g->SetLineWidth(0);
	// Double_t p0 = f1->GetParameter(0);
	// Double_t p1 = f1->GetParameter(1);
	// Double_t p2 = f1->GetParameter(2);
	// XNormPar[x*4] = p0;
	// XNormPar[x*4+1] = p1;
	// XNormPar[x*4+2] = p2;
	TF1 *f1 = new TF1("f1", "pol2", 200, 8000);
	f1->SetParameter(0, 0.0);
	f1->SetParameter(1, 1.0);
	f1->SetParameter(2, 0.0);
	g->Fit(f1, "q, rob=0.7+");
	g->SetMarkerStyle(7);
	g->SetLineWidth(0);
	Double_t p0 = f1->GetParameter(0);
	Double_t p1 = f1->GetParameter(1);
	Double_t p2 = f1->GetParameter(2);
	XNormPar[x*4] = p0;
	XNormPar[x*4+1] = p1;
	XNormPar[x*4+2] = p2;

	// residual
	TH2F *res = hresX[x];
	res->Reset();
	for (Int_t i = 0; i != pointNum; ++i) {
		Double_t px = g->GetPointX(i);
		Double_t py = g->GetPointY(i);
		res->Fill(f1->Eval(px)-py, py);
	}

	printf("%3s %3d %7lf %7lf %7lf\n", "X", x, p0, p1, p2);
	return 0;
}

int DSSD::NormYStrip(int xmin, int xmax, int y) {
	TGraph *g = gFitY[y];
	Int_t pointNum = 0;
	for (int i = xmin; i != xmax; ++i) {
		if (!xNormFlag[i]) continue;
		// // 选择x面需要的条
		// TCut xStripCut = TString::Format("xs==%d", i).Data();
		// // 选择y面需要的条
		// TCut yStripCut = TString::Format("ys==%d", y).Data();
		// // 需要画出的表达式，给参考条乘上刻度系数
		// TString varStr, hName;
		// hName.Form("h2Y[%d]X[%d]", y, i);
		// varStr.Form("%lf+%lf*xe+%lf*xe*xe:ye>>%s(1000, 0, 8000, 1000, 0, 8000)", XNormPar[i*4], XNormPar[i*4+1], XNormPar[i*4+2], hName.Data());
		// ipt->Draw(varStr.Data(), hit1or2Cut && xStripCut && yStripCut, "goff");

		// // add points to the graph
		// Long64_t rows = ipt->GetSelectedRows();
		// Double_t *pointY = ipt->GetV1();
		// Double_t *pointX = ipt->GetV2();
		// for (Long64_t j = 0; j != rows; ++j) {
		// 	g->SetPoint(pointNum, pointX[j], pointY[j]);
		// 	++pointNum;
		// }

		// if (LogFile) {
		// 	TH2F *h2xy = (TH2F*)gDirectory->Get(hName.Data());
		// 	LogFile->cd();
		// 	LogFile->GetDirectory("TH2")->cd();
		// 	h2xy->Write();
		// }

		TString gName;
		gName.Form("gX[%d]Y[%d]", i, y);
		if (!GraphFile) {
			printf("Error no GraphFile\n");
			return -1;
		}
		TGraph *gs = (TGraph*)GraphFile->Get((Name+"/"+gName).Data());
		Int_t n = gs->GetN();
		for (Int_t pi = 0; pi != n; ++pi) {
			Double_t px = gs->GetPointY(pi);
			Double_t py = gs->GetPointX(pi);
			py = XNormPar[i*4] + XNormPar[i*4+1]*py + XNormPar[i*4+2]*py*py;
			g->SetPoint(pointNum++, px, py);
		}
	}

	if (pointNum == 0) return -1;

	// fit
	TF1 *f1 = new TF1("f1", "pol2", 200, 8000);
	f1->SetParameter(0, 0.0);
	f1->SetParameter(1, 1.0);
	f1->SetParameter(2, 0.0);
	g->Fit(f1, "Q rob=0.7+");
	g->SetMarkerStyle(7);
	g->SetLineWidth(0);
	Double_t p0 = f1->GetParameter(0);
	Double_t p1 = f1->GetParameter(1);
	Double_t p2 = f1->GetParameter(2);
	YNormPar[y*4] = p0;
	YNormPar[y*4+1] = p1;
	YNormPar[y*4+2] = p2;


	// residual
	TH2F *res = hresY[y];
	res->Reset();
	for (Int_t i = 0; i != pointNum; ++i) {
		Double_t px = g->GetPointX(i);
		Double_t py = g->GetPointY(i);
		res->Fill(f1->Eval(px)-py, py);
	}


	printf("%3s %3d %7lf %7lf %7lf\n", "Y", y, p0, p1, p2);
	return 0;
}



MDSSD::MDSSD(TTree *iipt, TTree *oopt) {
	dssdNum = 0;
	opt = oopt;
	ipt = iipt;
	LogFile = nullptr;
}


MDSSD::~MDSSD() {
}


int MDSSD::Norm() {
	for (int i = 0; i != dssdNum; ++i) {
		dssds[i]->NormAll();
	}
	fill();
	return 0;
}


int MDSSD::AddDSSD(const char *name, int xs, int ys) {
	dssds.push_back(new DSSD(ipt, name, xs, ys));
	dssds[dssdNum]->LogFile = LogFile;
	LogFile->cd();
	TDirectoryFile *subDir = new TDirectoryFile(dssds[dssdNum]->Name.data(), dssds[dssdNum]->Name.data());
	subDir->cd();
	new TDirectoryFile("Fit", "Fit");
	new TDirectoryFile("Residual", "Residual");

	dssds[dssdNum]->GraphFile = GraphFile;
	++dssdNum;
	return 0;
}

int MDSSD::AddNode(int index, int xmin, int xmax, int ymin, int ymax, int dir) {
	dssds[index]->AddNode(xmin, xmax, ymin, ymax, dir);
	return 0;
}

int MDSSD::fill() {
	// calculate the sum of strip
	// int *sumXStrip = new int[dssdNum+1];
	// int *sumYStrip = new int[dssdNum+1];
	// sumXStrip[0] = sumYStrip[0] = 0;
	// for (int i = 1; i != dssdNum+1; ++i) {
	// 	sumXStrip[i] = sumXStrip[i-1] + dssds[i-1]->XStrip;
	// 	sumYStrip[i] = sumYStrip[i-1] + dssds[i-1]->YStrip;
	// }
	// // data to fill
	// Int_t *xhit = new Int_t[dssdNum];
	// Double_t *nxe = new Double_t[sumXStrip[dssdNum]];
	// Int_t *xe = new Int_t[sumXStrip[dssdNum]];
	// Int_t *xs = new Int_t[sumXStrip[dssdNum]];
	// Double_t *xres = new Double_t[sumXStrip[dssdNum]];

	// Int_t *yhit = new Int_t[dssdNum];
	// Double_t *nye = new Double_t[sumYStrip[dssdNum]];
	// Int_t *ye = new Int_t[sumYStrip[dssdNum]];
	// Int_t *ys = new Int_t[sumYStrip[dssdNum]];
	// Double_t *yres = new Double_t[sumYStrip[dssdNum]];
	vector<Int_t> xhit, yhit;
	vector<Double_t*> nxe, nye;
	vector<Int_t*> xe, ye, xs, ys;
	for (int i = 0; i != dssdNum; ++i) {
		xhit.push_back(0);
		yhit.push_back(0);
		xe.push_back(new Int_t[dssds[i]->XStrip]);
		nxe.push_back(new Double_t[dssds[i]->XStrip]);
		ye.push_back(new Int_t[dssds[i]->YStrip]);
		nye.push_back(new Double_t[dssds[i]->YStrip]);
		xs.push_back(new Int_t[dssds[i]->XStrip]);
		ys.push_back(new Int_t[dssds[i]->XStrip]);
	}

	// branch
	for (int i = 0; i != dssdNum; ++i) {
		DSSD* dssd = dssds[i];
		// name for branch
		string xhitStr = dssd->Name+"xhit";
		string yhitStr = dssd->Name+"yhit";
		string xeStr = dssd->Name+"xe";
		string yeStr = dssd->Name+"ye";
		string xsStr = dssd->Name+"xs";
		string ysStr = dssd->Name+"ys";

	// output tree branch
		opt->Branch(xhitStr.data(), &xhit[i], (xhitStr+"/I").data());
		opt->Branch(xeStr.data(), nxe[i], (xeStr+"["+xhitStr+"]/D").data());
		opt->Branch(xsStr.data(), xs[i], (xsStr+"["+xhitStr+"]/I").data());
		opt->Branch(yhitStr.data(), &yhit[i], (yhitStr+"/I").data());
		opt->Branch(yeStr.data(), nye[i], (yeStr+"["+yhitStr+"]/D").data());
		opt->Branch(ysStr.data(), ys[i], (ysStr+"["+yhitStr+"]/I").data());

		// input tree branch
		ipt->SetBranchAddress(xhitStr.data(), &xhit[i]);
		ipt->SetBranchAddress(xeStr.data(), xe[i]);
		ipt->SetBranchAddress(xsStr.data(), xs[i]);
		ipt->SetBranchAddress(yhitStr.data(), &yhit[i]);
		ipt->SetBranchAddress(yeStr.data(), ye[i]);
		ipt->SetBranchAddress(ysStr.data(), ys[i]);
	}

	// loop and fill
	printf("filling   0%%");
	fflush(stdout);
	Long64_t nentry = ipt->GetEntries();
	Long64_t nentry100 = nentry / 100;
	for (Long64_t jentry = 0; jentry != nentry; ++jentry) {
		ipt->GetEntry(jentry);
		for (int i = 0; i != dssdNum; ++i) {
			for (int j = 0; j != xhit[i]; ++j) {
				// int index = j + sumXStrip[i];
				Double_t p0 = dssds[i]->XNormPar[xs[i][j]*4];
				Double_t p1 = dssds[i]->XNormPar[xs[i][j]*4+1];
				Double_t p2 = dssds[i]->XNormPar[xs[i][j]*4+2];
				nxe[i][j] = p0 + p1*xe[i][j] + p2*xe[i][j]*xe[i][j];
				// nxe[i][j] = p0 + p1 * xe[i][j];
			}
			for (int j = 0; j != yhit[i]; ++j) {
				// int index = j + sumXStrip[i];
				Double_t p0 = dssds[i]->YNormPar[ys[i][j]*4];
				Double_t p1 = dssds[i]->YNormPar[ys[i][j]*4+1];
				Double_t p2 = dssds[i]->YNormPar[ys[i][j]*4+2];
				nye[i][j] = p0 + p1*ye[i][j] + p2*ye[i][j]*ye[i][j];
				// nye[i][j] = p0 + p1 * ye[i][j];
			}
		}
		opt->Fill();

		if (jentry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", jentry / nentry100);
			fflush(stdout);
		}
	}
	printf("\b\b\b\b100%%\n");

	for (int i = 0; i != dssdNum; ++i) {
		delete[] nxe[i];
		delete[] xe[i];
		delete[] xs[i];
		delete[] nye[i];
		delete[] ye[i];
		delete[] ys[i];
	}

	return 0;
}

int MDSSD::SetLogFile(TFile *f) {
	if (f == nullptr) return -1;
	LogFile = f;
	for (int i = 0; i != dssdNum; ++i) {
		dssds[i]->LogFile = f;
		f->cd();
		TDirectoryFile *subDir = new TDirectoryFile(dssds[i]->Name.data(), dssds[i]->Name.data());
		subDir->cd();
		new TDirectoryFile("Fit", "Fit");
		new TDirectoryFile("Residual", "Residual");
	}
	return 0;
}


int MDSSD::SetGraphFile(TFile *f) {
	if (f == nullptr) return -1;
	GraphFile = f;
	for (int i = 0; i != dssdNum; ++i) {
		dssds[i]->GraphFile = f;
	}
	return 0;
}
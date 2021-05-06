#include "Simulation.h"


Simulation::Simulation() {
	gr = new TRandom3(0);

	detectors = 0;
	particles = 0;
	targetSigma = nullptr;
	z = nullptr;
	batches = 0;

	logFile = nullptr;
	logTree = nullptr;
	minLoss = 100.0;
	minRes = nullptr;

	nodeN = 0;
	nodes.clear();

	results = nullptr;

	testNum = 0;
}

Simulation::~Simulation() {
	if (gr) delete gr;

	if (z) delete[] z;
	if (targetSigma) delete[] targetSigma;

	if (minRes) delete[] minRes;
	if (logFile) logFile->Close();

	if (batches) {
		for (int i = 0; i != batches; ++i) {
			delete[] results[i].Res;
		}
		delete[] results;
	}
}

int Simulation::SetDetectors(int detN, const Double_t *zz, const Double_t *ss) {
	detectors = detN;
	targetSigma = new Double_t[detectors];
	z = new Double_t[detectors];
	for (int j = 0; j != detectors; ++j) {
		z[j] = zz[j];
		targetSigma[j] = ss[j];
	}
	return 0;
}


int Simulation::SetTraceGen(TTree *ipt, const char* sk, const char *sb) {
	tree = ipt;
	tree->SetBranchAddress(sk, &rk);
	tree->SetBranchAddress(sb, &rb);
	return 0;
}


int Simulation::AddNode(int pp, int maxIter, Double_t eps, Double_t aa, Double_t rr) {
	if (pp > particles) particles = pp;
	SimNode node;
	node.Particles = pp;
	node.MaxIter = maxIter;
	node.Eps = eps;
	node.AlphaInit = aa;
	node.ReduceAlpha = rr;
	nodes.push_back(node);
	nodeN++;
	return 0;
}


int Simulation::Train(int bb) {
	batches = bb;
	// result
	results = new SimResult[batches];
	for (int i = 0; i != batches; ++i) {
		results[i].Loss = 100.0;
		results[i].Res = new Double_t[detectors];
		results[i].Particles = particles;
		results[i].Detectors = detectors;
	}

	// variables
	Double_t *rp = new Double_t[particles*detectors];	// “真实”的位置,rp[particle][detector]
	Double_t *res = new Double_t[detectors];			// 当前的探测器分辨
	TH1D **hdp = new TH1D*[detectors];					// 当前的残差分布
	Double_t *ns = new Double_t[detectors];				// 当前的残差分布的sigma


	// log
	TGraph *gLoss = nullptr;
	TGraph **gRes = nullptr;


	for (int batch = 0; batch != batches; ++batch) {
		printf("\n============================================================\n");
		printf("batch %d:\n", batch);

		SimResult &result = results[batch];
		if (logFile) {
			logFile->cd();
			TString dirName;
			dirName.Form("batch%03d", batch);
			TDirectoryFile *dirBatch = new TDirectoryFile(dirName.Data(), dirName.Data());
			dirBatch->cd();
			gLoss = new TGraph;						// loss 变化
			gLoss->SetName(TString::Format("loss%03d", batch));
			gRes = new TGraph*[detectors];			// res 变化
			for (int j = 0; j != detectors; ++j) {
				gRes[j] = new TGraph;
				gRes[j]->SetName(TString::Format("res%d_%03d", j, batch));
			}
		}


		// generate the trace and the position on the detectors
		genTrace(rp);
		// random init resolution
		for (int j = 0; j != detectors; ++j) {
			res[j] = gr->Uniform(0.2, 0.6);
			result.Res[j] = res[j];
		}


		currentNode = 0;
		for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
			alpha = iter->AlphaInit;
			for (int iteration = 0; iteration != iter->MaxIter; ++iteration) {
				hitAndTrace(rp, res, hdp);
				fitResidual(hdp, ns);
				Double_t loss = lossFunc(ns, iter->ReduceAlpha);
				if (logFile) {
					gLoss->SetPoint(iteration, iteration, loss);
					for (int j = 0; j != detectors; ++j) {
						gRes[j]->SetPoint(iteration, iteration, res[j]);
					}
				}
				if (loss < result.Loss) {
					result.Loss = loss;
					for (int j = 0; j != detectors; ++j) {
						result.Res[j] = res[j];
					}
				}
				if (loss < iter->Eps) break;
				step(ns, res);
				printf("batch %d, node %d, iteration %d, loss %lf\n", batch, currentNode, iteration, loss);
			}
			currentNode++;
			printf("------------------------------------------------------------\n");
		}


		// print result
		printf("\n");
		for (int j = 0; j != detectors-1; ++j) {
			printf("%lf, ", result.Res[j]);
		}
		printf("%lf\n\n", result.Res[detectors-1]);


		if (logFile) {
			if (gLoss) {
				gLoss->Write();
				delete gLoss;
			}
			if (gRes) {
				for (int j = 0; j != detectors; ++j) {
					gRes[j]->Write();
					delete gRes[j];
				}
				delete[] gRes;
			}
			if (logTree) {
				minLoss = result.Loss;
				for (int j = 0; j != detectors; ++j) {
					minRes[j] = result.Res[j];
				}
				logTree->Fill();
			}
		}
	}

	if (logTree) {
		logFile->cd();
		logTree->Write();
	}

	delete[] res;
	delete[] hdp;
	delete[] rp;
	delete[] ns;

	for (int i = 0; i != batches; ++i) {
		printf("\nLoss %lf,  Particles %d,  Detecotrs %d\n", results[i].Loss, results[i].Particles, results[i].Detectors);
		printf("Res\n");
		for (int j = 0; j != detectors; ++j) {
			printf("%lf ", results[i].Res[j]);
		}
		printf("\b\n");
	}

	return 0;
}




int Simulation::Test(int bn, int cnt) {
	Double_t *rp = new Double_t[particles*detectors];	// “真实”的位置,rp[particle][detector]
	TH1D **hdp = new TH1D*[detectors];					// 当前的残差分布
	Double_t *ns = new Double_t[detectors];				// 当前的残差分布的sigma
	Double_t loss;

	TTree *opt = nullptr;
	TDirectoryFile *testDir = nullptr;
	if (logFile) {
		TString dirName;
		dirName.Form("test%d", testNum);
		testDir = new TDirectoryFile(dirName.Data(), dirName.Data());
		testDir->cd();
		opt = new TTree("tree", "log of test");
		opt->Branch("loss", &loss, "loss/D");
		opt->Branch("sigma", ns, TString::Format("sigma[%d]/D", detectors));
	}

	AddNode(particles, 0, 1e-6, 0.0, 0.0);
	for (int i = 0; i != cnt; ++i) {
		printf("\n============================================================\n");
		printf("batch %d, test %d:\n", bn, i);
		if (logFile) {
			logFile->cd();
			if (testDir) testDir->cd();
			if (opt) opt->Fill();
			TString dirName;
			dirName.Form("subTest%02d", i);
			TDirectoryFile *subDir = new TDirectoryFile(dirName.Data(), dirName.Data());
			subDir->cd();
		}

		genTrace(rp);
		hitAndTrace(rp, results[bn].Res, hdp);
		// save hdp
		if (logFile) {
			for (int j = 0; j != detectors; ++j) {
				hdp[j]->Write();
			}
		}
		fitResidual(hdp, ns);
		loss = lossFunc(ns);

		printf("\ntest %d loss %lf\n\n", i, loss);
		printf("| target   | present  | resolution |\n");
		printf("| -------- | -------- | -----------|\n");
		for (int j = 0; j != detectors; ++j) {
			printf("| %lf | %lf | %lf   |\n", targetSigma[j], ns[j], results[bn].Res[j]);
		}

	}

	if (logFile && testDir && opt) {
		logFile->cd();
		testDir->cd();
		opt->Write();
	}



	delete[] hdp;
	delete[] ns;
	delete[] rp;

	testNum++;
	return 0;
}

int Simulation::SetLog(const char *opf) {
	if (!opf) return -1;
	if (!detectors) return -1;
	logFile = new TFile(opf, "recreate");
	logTree = new TTree("tree", "log of training");
	minRes = new Double_t[detectors];
	logTree->Branch("res", minRes, TString::Format("res[%d]/D", detectors));
	logTree->Branch("loss", &minLoss, "loss/D");
	return 0;
}


// 模拟生成径迹并计算径迹和探测器的交点
int Simulation::genTrace(Double_t *rp) {
	// log trace 记录“真实”径迹
	TH2I *hTrace = nullptr;
	if (logFile) {
		hTrace = new TH2I("hTrace", "trace for 'real' particles", 220, -2000.0, 200.0, 300, -150, 150);
	}

	Long64_t nentries = tree->GetEntriesFast();
	Long64_t jentry = nentries * gr->Uniform(0, 1);
	fprintf(stderr, "generating    0%%");
	for (int i = 0; i != particles; ++i) {
		tree->GetEntry(jentry);
		rk = gr->Gaus(rk, rk/300.0);
		rb = gr->Gaus(rb, rb/300.0);
		for (int j = 0; j != detectors; ++j) {
			rp[i*detectors+j] = rk * z[j] + rb;			// “真实”位置
		}
		jentry++;
		jentry = jentry > nentries ? 0 : jentry;
		if (logFile) {
			addTrace(hTrace, rk, rb, -1800.0, 0.0);
		}
		if ((i+1) % (particles/100) == 0) fprintf(stderr, "\b\b\b\b%3d%%", (i+1) / (particles/100));
		fflush(stderr);
	}
	fprintf(stderr, "\n");
	if (logFile) {
		hTrace->Write();
		if (hTrace) delete hTrace;
	}
	return 0;
}


// 线性拟合
Double_t Simulation::sFit(Double_t *fx, Double_t *fy, Double_t &k, Double_t &b) {
	int n = detectors;
	Double_t sumx = 0.0;
	Double_t sumy = 0.0;
	Double_t sumxy = 0.0;
	Double_t sumx2 = 0.0;
	for (int i = 0; i != n; ++i) {
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
	for (int i = 0; i != n; ++i) {
		Double_t x = fx[i];
		Double_t y = fy[i];
		Double_t t = y - k*x - b;
		chi2 += t * t;
	}
	return chi2;
}


// 模拟粒子打在探测器上，根据分辨产生信号，根据信号重建径迹
int Simulation::hitAndTrace(Double_t *rp, Double_t *res, TH1D **hdp) {
	for (int j = 0; j != detectors; ++j) {
		TString name;
		name.Form("hdp%d", j);
		hdp[j] = new TH1D(name.Data(), name.Data(), 200, -1, 1);
	}

	int pp = nodes[currentNode].Particles;
	// hit and trace
	Double_t *vp = new Double_t[detectors];
	for (int i = 0; i != pp; ++i) {
		for (int j = 0; j != detectors; ++j) {
			vp[j] = gr->Gaus(rp[i*detectors+j], res[j]);			// hit and bring resolution
		}
		Double_t k, b;
		sFit(z, vp, k, b);
		// if (sFit((Double_t*)z, (Double_t*)vp, k, b) >= CHI_THRESHOLD) continue;			// no need, the chi2 is very small
		for (int j = 0; j != detectors; ++j) {
			Double_t tp = b + k * z[j];			// 用径迹拟合的位置
			hdp[j]->Fill(vp[j]-tp);					// 填充残差
		}
	}
	delete[] vp;
	return 0;
}



// 拟合残差，hdp被拟合的分布，ns拟合结果sigma
int Simulation::fitResidual(TH1D **hdp, Double_t *ns) {
	TF1 *f1 = new TF1("f1", "gaus", -0.6, 0.6);
	for (int j = 0; j != detectors; ++j) {
		hdp[j]->Fit(f1, "QRN+");
		ns[j] = f1->GetParameter(2);
	}
	delete f1;
	for (int j = 0; j != detectors; ++j) {
		delete hdp[j];
	}
	return 0;
}



// 计算loss函数，ns - 当前残差sigma，rr - 学习率衰变因子
Double_t Simulation::lossFunc(Double_t *ns, Double_t rr) {
	Double_t l = 0.0;
	for (int j = 0; j != detectors; ++j) {
		Double_t t = ns[j] - targetSigma[j];
		l += t * t;
	}
	if (l > lastLoss) {
		alpha *= rr;
	}
	lastLoss = l;
	return l;
}


// 根据结果修改分辨率
int Simulation::step(Double_t *ns, Double_t *res) {
	for (int j = 0; j != detectors; ++j) {
		Double_t t = targetSigma[j] - ns[j];
		res[j] += alpha * t;
	}
	return 0;
}


int Simulation::addTrace(TH2 *h, Double_t k, Double_t b, Double_t minBin, Double_t maxBin) {
	if (!h) return -1;
	if (minBin >= maxBin) return -1;
	const int step = 10;
	for (int i = minBin; i < maxBin; i += step) {
		h->Fill(i, Int_t(i*k+b));
	}
	return 0;
}

#include "Simulation.h"
#include <ctime>

// detector data
const int DETECTORS = 5;
const int PARTICLES = 10000000;
// position
const Double_t xz[DETECTORS]{-1750.7, -1713.3, -1250.7, -1213.3, -381.6};
const Double_t yz[DETECTORS]{-1742.1, -1721.9, -1242.1, -1221.9, -373.4};
// resolution
const Double_t sx[DETECTORS]{0.336, 0.365, 0.399, 0.405, 0.2027};
const Double_t sy[DETECTORS]{0.365, 0.375, 0.425, 0.415, 0.2349};
// step
const int MAXRUN = 100;
const double eps = 1e-6;
const int MAXBATCH = 10;


int main(int argc, char **argv) {
	TFile *ipf = new TFile("../../data/tracking.root");
	TTree *ipt = (TTree*)ipf->Get("tree");

	Simulation *xSim = new Simulation;
	xSim->SetDetectors(DETECTORS, xz, sx);
	xSim->SetTraceGen(ipt, "kx", "bx");
	xSim->SetLog("../../log/xSimLog.root");
	xSim->AddNode(PARTICLES/100, 20, 1e-4, 1.0, 0.8);
	xSim->AddNode(PARTICLES/10, 40, 1e-5, 0.8, 0.8);
	xSim->AddNode(PARTICLES, 10, 1e-6, 0.1, 0.8);
	clock_t t = clock();
	xSim->Train(10);
	printf("%f s\n", float(clock()-t)/CLOCKS_PER_SEC);
	// xSim->Test(0, 3);
	delete xSim;

	Simulation *ySim = new Simulation;
	ySim->SetDetectors(DETECTORS, yz, sy);
	ySim->SetTraceGen(ipt, "ky", "by");
	ySim->SetLog("../../log/ySimLog.root");
	ySim->AddNode(PARTICLES/100, 20, 1e-4, 1.0, 0.8);
	ySim->AddNode(PARTICLES/10, 40, 1e-5, 0.8, 0.8);
	ySim->AddNode(PARTICLES, 10, 1e-6, 0.1, 0.8);
	ySim->Train(10);
	// ySim->Test(10);
	delete ySim;

	ipf->Close();
	return 0;
}
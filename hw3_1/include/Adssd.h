#ifndef __ADSSD_H__
#define __ADSSD_H__

#include <vector>
#include <queue>
#include <map>
#include <functional>
#include <algorithm>

#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include <TLatex.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TFitResult.h>
#include <TArrayD.h>

#include "AdssdBase.h"

class Adssd: public AdssdBase {
public:
	Adssd(TTree *tree, TFile *f = nullptr);
	~Adssd();


	virtual int Analysis(TTree *oopt = nullptr);
	virtual int SetAlpha(int an, const Double_t *ap);
private:
	const int chp;

	TRandom3 *gr;

	Int_t alphaPeakNum;
	const Double_t *alphaPeaks;

	TSpectrum *spectrum;

	TFile *logFile;

	TTree *opt;

	Double_t *ccpe;

	virtual int selectPeaks(std::vector<Double_t> &spe, std::vector<Double_t> &pe);
	virtual int peaks(TH1F *h, std::vector<Double_t> &pe, Double_t thres = 0.05, int backsub = 0);
	virtual int changePeaks(int *cpe, int peaks);
	virtual int gaussFit(TH1F *h, const std::vector<Double_t> &spe, Double_t *par);
	virtual int calibration(Double_t *fitPar, int ch, Double_t &p0, Double_t &p1, Double_t &chi2);
	virtual int setBranch();
};



#endif
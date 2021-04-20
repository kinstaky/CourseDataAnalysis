#ifndef Tracking_h
#define Tracking_h

#include "trackingBase.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH2D.h>
#include <TFitResult.h>

class Tracking: public trackingBase {
public :
	Tracking(TTree *tree = 0);
	virtual ~Tracking();
	virtual void Loop();
private:
    // anode
    Double_t anode[5];
    // variable for x tracking
    Double_t xx[5], xz[5];
    // variable for y tracking
    Double_t yy[5], yz[5];
    // residual
    Double_t dx[5], dy[5];
    // target position
    Double_t tx, ty, tz;
    // target position after projection
    Double_t ptx, pty;
    // chi2/ndf
    Double_t c2nx, c2ny;
    // number of PPAC to track, 0 for tracking failure
    // 用于重建的探测器数量，0表示不满足重建条件
    Int_t numX, numY;
    // combination of the valid PPAC,用于重建径迹的探测器组合
    // 使用标志位的方式记录
    Int_t trackX, trackY;
    // 拟合得到的径迹方程
    Double_t bx, kx;
    Double_t by, ky;
    // 重建后的径迹在不同探测器上的位置
    Double_t fxx[5], fyy[5];

    // init the variables, return: 0 -- can't track, 1 -- track x, 2 -- track y, 3 -- track x and y
    virtual Int_t trackInit();
    virtual void setBranch(TTree *tree);
    virtual void addTrace(TH2D *h, Double_t k, Double_t b, Int_t minBin, Int_t maxBin);
};

#endif
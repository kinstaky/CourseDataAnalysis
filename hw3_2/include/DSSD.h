#ifndef __DSSD_H__
#define __DSSD_H__

#include <TTree.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCut.h>
#include <TH2F.h>

#include <vector>
#include <string>


const int DIR_X2Y = 0;
const int DIR_Y2X = 1;


struct NodeInfo {
	int direction;			// 0: from x to y 	1: from y to x
	int xmin;
	int xmax;
	int ymin;
	int ymax;
};


class DSSD {
public:
	DSSD(TTree *iipt, const char *nn, int xs, int ys);			// nn - name, xs - xstrip, ys - ystrip
	~DSSD();

	// 添加节点，即刻度的过程。规定了怎么刻度，即用哪一面的哪几根条刻度另一面的哪几根条
	int AddNode(int xmin, int xmax, int ymin, int ymax, int dir);
	// 刻度所有条
	int NormAll();
	// 刻度x面
	int NormXStrip(int ymin, int ymax, int x);
	// 刻度y面
	int NormYStrip(int xmin, int xmax, int y);

	// normalize parameter, [i*2] -- k of strip[i], [i*2+1] -- b for strip[i]
	// normalize energy = k * energy + b
	// 刻度系数
	std::vector<Double_t> XNormPar;
	std::vector<Double_t> YNormPar;
	// dssd 的名字，用于在TTree中读取变量
	std::string Name;
	// x 和 y 方向条的总数
	int XStrip, YStrip;

	TFile *LogFile;
	TFile *GraphFile;
private:
	TTree *ipt;

	// nodes即刻度的节点，一个节点代表一个步骤
	std::vector<NodeInfo> nodes;
	// TCut，用于选择单举和双举事件
	// TCut hit1Cut, xhit2Cut, yhit2Cut, hit2Cut, hit1or2Cut;
	// 归一化标志位
	std::vector<bool> xNormFlag;
	std::vector<bool> yNormFlag;
	// 归一化残差
	std::vector<TH2F*> hresX;
	std::vector<TH2F*> hresY;
	// 归一化拟合
	std::vector<TGraph*> gFitX;
	std::vector<TGraph*> gFitY;
	// 保存刻度系数的TGraph
	TGraph *gxNormP[3];
	TGraph *gyNormP[3];

	int genGraph();
};


class MDSSD {
public:
	MDSSD(TTree *ipt, TTree* oopt);
	~MDSSD();

	int Norm();
	int AddDSSD(const char *name, int xs, int ys);
	int AddNode(int index, int xmin, int xmax, int ymin, int ymax, int dir);
	int SetLogFile(TFile *f);
	int SetGraphFile(TFile *f);

	TFile *LogFile;
	TFile *GraphFile;
private:
	int dssdNum;
	std::vector<DSSD*> dssds;
	TTree *opt;
	TTree *ipt;

	int fill();
};


#endif
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <vector>
#include <TH2I.h>

using std::vector;

// 模拟结果类
struct SimResult {
	Double_t Loss;
	Double_t *Res;
	int Particles;
	int Detectors;
};


// 模拟模型节点数据
struct SimNode {
	// physics
	int Particles;
	// iteration
	int MaxIter;
	Double_t Eps;
	// learning
	Double_t AlphaInit;
	Double_t ReduceAlpha;
};



// 模拟器类
class Simulation {
public:
	Simulation();
	~Simulation();
	// 设置探测器,
	// 	detN	探测器数量
	//	zz		探测器位置
	//	ss 		探测器目标残差
	int SetDetectors(int detN, const Double_t *zz, const Double_t *ss);
	// 设置径迹生成参数
	//	ipt 	径迹的树
	// 	sk 		树中径迹斜率k的branch名字
	//	sb 		树中径迹截距b的branch名字
	int SetTraceGen(TTree *ipt, const char* sk, const char *sb);
	// 设置输出
	//	fileName 输出文件名
	int SetLog(const char *fName = nullptr);
	// 设置模型节点
	// 	pp		粒子数量
	// 	maxIter	迭代次数
	// 	eps		迭代精度
	// 	aa 		初始学习率alpha
	// 	rr 		学习率衰减因子
	int AddNode(int pp, int maxIter, Double_t eps, Double_t aa, Double_t rr);
	// 开始拟合训练
	//  bb		训练组数
	int Train(int bb = 1);
	// 测试结果
	// 	bn		训练结果序号
	//  opf		输出文件
	//  cnt 	测试次数
	int Test(int bn, int cnt = 1);
private:
	int detectors;
	int particles;
	int batches;

	TRandom3 *gr;					// 随机数引擎
	TTree *tree;					// 用于随机径迹的数据
	Double_t rk, rb;				// “真实”径迹的k和b
	Double_t *z;					// 探测器的z平面
	Double_t *targetSigma;			// 目标残差的sigma

	vector<SimNode> nodes;			// 模型节点
	int nodeN;						// 节点总数
	int currentNode;				// 当前使用节点

	Double_t lastLoss;
	Double_t alpha;

	TFile *logFile;					// 记录的文件
	TTree *logTree;					// 记录的树
	Double_t minLoss;				// 记录的最小平方差
	Double_t *minRes;				// 记录的最接近的分辨

	int testNum;					// test id

	SimResult* results;				// 结果

	// 生成径迹
	int genTrace(Double_t *rp);
	// simple linear fit, return the chi2,
	Double_t sFit(Double_t *fx, Double_t *fy, Double_t &k, Double_t &b);
	// 模拟粒子打在探测器上，因为分辨位置偏移，利用偏移后的位置重建径迹，返回残差分布
	int hitAndTrace(Double_t *rp, Double_t *res, TH1D **hdp);
	// 拟合残差
	int fitResidual(TH1D **hdp, Double_t *ns);
	// 计算当前迭代的残差和目标残差的平方差
	Double_t lossFunc(Double_t *ns, Double_t rr = 0.8);
	// 根据学习率修改分辨用以下次迭代
	int step(Double_t *ns, Double_t *res);
	// 记录径迹
	int addTrace(TH2 *h, Double_t k, Double_t b, Double_t minBin = -1800.0, Double_t maxBin = 0.0);
};
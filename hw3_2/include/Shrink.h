#ifndef __SHRINK_H__
#define __SHRINK_H__

#include <vector>
#include <map>

#include <TTree.h>
#include <TString.h>


class DSSDShrink {
public:
	DSSDShrink(TTree *it, TTree *st, TTree *rt);
	~DSSDShrink();

	int AddDSSD(const char *nn, int xs, int ys, int thre = 100);
	int Shrink();
private:
	std::vector<TString> names;
	std::vector<int> xStrip;
	std::vector<int> yStrip;
	std::vector<int> threshold;
	int dssdNum;

	TTree *ipt, *stree, *rtree;

	struct Event{
		Double_t Xe;
		Double_t Ye;
		Double_t XStrip;
		Double_t YStrip;
		Int_t Flag;
	};
};

#endif
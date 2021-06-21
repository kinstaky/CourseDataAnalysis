#ifndef __DSSD_H__
#define __DSSD_H__

#include <map>

#include <TTree.h>


class DSSD {
public:
	DSSD(TTree *iipt, TTree *oopt, Int_t fs, Int_t bs);
	~DSSD();

	int Analysis();

private:
	TTree *ipt, *opt;
	Int_t fStrip, bStrip;


	struct dssdInfo {
		Double_t Energy;
		Int_t Strip;
		Int_t Hit;
		Int_t Index;				// hit index
	};

	struct dsfbInfo {
		ULong64_t Time;
		Double_t FStrip, BStrip;
		Double_t FEnergy, BEnergy;
	};

	std::multimap<ULong64_t, dssdInfo> mapfs, mapbs;
	std::multimap<ULong64_t, dsfbInfo> mapfbs;

	Double_t re;			// raw energy
	ULong64_t rt;			// raw time
	Int_t rs, side;		// raw strip, side

	Int_t chit;
	Double_t cfStrip[32], cbStrip[32];
	Double_t cEnergy[32], cfe[32], cbe[32];
	ULong64_t cTime[32];

	int setBranch();
	int fillMap();
	int printMap(int entry = 20);
	int searchWindow(ULong64_t tw = 100000, ULong64_t toff = 0);
	int hecorr(int tw, int toff);
	int getHit(int tw, int toff);
	int correlation(int tw, int toff);
};

#endif
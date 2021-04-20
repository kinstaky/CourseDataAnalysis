//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 19 14:52:03 2021 by ROOT version 6.23/01
// from TTree tree/tree
// found on file: ../data/f8ppac001.root
//////////////////////////////////////////////////////////

#ifndef trackingBase_h
#define trackingBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class trackingBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         PPACF8[5][5];
   Float_t         F8PPACRawData[5][5];
   Int_t           beamTrig;
   Int_t           must2Trig;
   Float_t         targetX;
   Float_t         targetY;

   // List of branches
   TBranch        *b_PPACF8;   //!
   TBranch        *b_F8PPACRawData;   //!
   TBranch        *b_beamTrig;   //!
   TBranch        *b_must2Trig;   //!
   TBranch        *b_targetX;   //!
   TBranch        *b_targetY;   //!

   trackingBase(TTree *tree=0);
   virtual ~trackingBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef trackingBase_cxx
trackingBase::trackingBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/f8ppac001.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/f8ppac001.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

trackingBase::~trackingBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trackingBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trackingBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void trackingBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PPACF8", PPACF8, &b_PPACF8);
   fChain->SetBranchAddress("F8PPACRawData", F8PPACRawData, &b_F8PPACRawData);
   fChain->SetBranchAddress("beamTrig", &beamTrig, &b_beamTrig);
   fChain->SetBranchAddress("must2Trig", &must2Trig, &b_must2Trig);
   fChain->SetBranchAddress("targetX", &targetX, &b_targetX);
   fChain->SetBranchAddress("targetY", &targetY, &b_targetY);
   Notify();
}

Bool_t trackingBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trackingBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trackingBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trackingBase_cxx

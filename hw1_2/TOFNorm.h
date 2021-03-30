//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 30 20:59:45 2021 by ROOT version 6.23/01
// from TTree tree/ADC tree with veto wall
// found on file: AddVeto.root
//////////////////////////////////////////////////////////

#ifndef TOFNorm_h
#define TOFNorm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class TOFNorm {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        x;
   Double_t        e;
   Double_t        tof;
   Int_t           pid;
   Double_t        tu;
   Double_t        td;
   Double_t        qu;
   Double_t        qd;
   Int_t           itu;
   Int_t           itd;
   Int_t           iqu;
   Int_t           iqd;
   Double_t        diff;
   Double_t        vx;
   Double_t        vtof;
   Double_t        vtu;
   Double_t        vtd;
   Double_t        vqu;
   Double_t        vqd;
   Int_t           vitu;
   Int_t           vitd;
   Int_t           viqu;
   Int_t           viqd;

   // List of branches
   TBranch        *b_x;   //!
   TBranch        *b_e;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_tu;   //!
   TBranch        *b_td;   //!
   TBranch        *b_qu;   //!
   TBranch        *b_qd;   //!
   TBranch        *b_itu;   //!
   TBranch        *b_itd;   //!
   TBranch        *b_iqu;   //!
   TBranch        *b_iqd;   //!
   TBranch        *b_diff;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vtof;   //!
   TBranch        *b_vtu;   //!
   TBranch        *b_vtd;   //!
   TBranch        *b_vqu;   //!
   TBranch        *b_vqd;   //!
   TBranch        *b_vitu;   //!
   TBranch        *b_vitd;   //!
   TBranch        *b_viqu;   //!
   TBranch        *b_viqd;   //!

   TOFNorm(TTree *tree=0);
   virtual ~TOFNorm();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TOFNorm_cxx
TOFNorm::TOFNorm(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AddVeto.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("AddVeto.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

TOFNorm::~TOFNorm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TOFNorm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TOFNorm::LoadTree(Long64_t entry)
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

void TOFNorm::Init(TTree *tree)
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

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("e", &e, &b_e);
   fChain->SetBranchAddress("tof", &tof, &b_tof);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("tu", &tu, &b_tu);
   fChain->SetBranchAddress("td", &td, &b_td);
   fChain->SetBranchAddress("qu", &qu, &b_qu);
   fChain->SetBranchAddress("qd", &qd, &b_qd);
   fChain->SetBranchAddress("itu", &itu, &b_itu);
   fChain->SetBranchAddress("itd", &itd, &b_itd);
   fChain->SetBranchAddress("iqu", &iqu, &b_iqu);
   fChain->SetBranchAddress("iqd", &iqd, &b_iqd);
   fChain->SetBranchAddress("diff", &diff, &b_diff);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vtof", &vtof, &b_vtof);
   fChain->SetBranchAddress("vtu", &vtu, &b_vtu);
   fChain->SetBranchAddress("vtd", &vtd, &b_vtd);
   fChain->SetBranchAddress("vqu", &vqu, &b_vqu);
   fChain->SetBranchAddress("vqd", &vqd, &b_vqd);
   fChain->SetBranchAddress("vitu", &vitu, &b_vitu);
   fChain->SetBranchAddress("vitd", &vitd, &b_vitd);
   fChain->SetBranchAddress("viqu", &viqu, &b_viqu);
   fChain->SetBranchAddress("viqd", &viqd, &b_viqd);
   Notify();
}

Bool_t TOFNorm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TOFNorm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TOFNorm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TOFNorm_cxx

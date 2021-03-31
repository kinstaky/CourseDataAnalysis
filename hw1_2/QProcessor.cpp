#define QProcessor_cxx
#include "QProcessor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void QProcessor::Loop()
{
//   In a ROOT session, you can do:
//      root> .L QProcessor.C
//      root> QProcessor t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;

    // new file
    TFile *opf = new TFile("QProcessor.root", "recreate");
    TTree *opt = new TTree("tree", "q process step 1");
    // new data
    Double_t q0;
    // new branch
    opt->Branch("q0", &q0, "q0/D");

    fChain->AddFriend(opt);

    // qu
    fChain->Draw("iqu >> hqu(420, 0, 4200)", "", "goff");
    TH1D *hqu = (TH1D*)gDirectory->Get("hqu");
    TF1 *f1 = new TF1("f1", "gaus", 0, 250);
    hqu->Fit(f1, "R");
    Double_t pedUSigma = f1->GetParameter(2);
    Double_t pedU = f1->GetParameter(1);
    TString sPedU(Form("%lf", pedU));
    TString squLimit(Form("iqu > %lf+3.0*%lf && iqu < 4095", pedU, pedUSigma));
    opt->SetAlias("qua", ("iqu-"+sPedU).Data());
    opt->SetAlias("quLimit", squLimit.Data());
    delete f1;

    // qd
    fChain->Draw("iqd >> hqd(420, 0, 4200)", "", "goff");
    TH1D *hqd = (TH1D*)gDirectory->Get("hqd");
    f1 = new TF1("f1", "gaus", 0, 250);
    hqd->Fit(f1, "R");
    Double_t pedDSigma = f1->GetParameter(2);
    Double_t pedD = f1->GetParameter(1);
    TString sPedD{Form("%lf", pedD)};
    TString sqdLimit(Form("iqd > %lf+3.0*%lf && iqd < 4095", pedD, pedDSigma));
    opt->SetAlias("qda", ("iqd-"+sPedD).Data());
    opt->SetAlias("qdLimit", sqdLimit.Data());
    opt->SetAlias("qLimit", "qdLimit && quLimit");
    delete f1;

    // fChain->Draw("qda:qua>>(200,0,4000,200,0,4000)", "!qLimit", "colz");

    Long64_t nentries = fChain->GetEntries();
    for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
        fChain->GetEntry(jentry);
        q0 = sqrt((iqu-pedUSigma) * (iqd-pedDSigma));
        opt->Fill();
    }

    fChain->Draw("q0:(itd+itu)/2.0", "qLimit", "colz");

    opt->Write();
    opf->Close();
}

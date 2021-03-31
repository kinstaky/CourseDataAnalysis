#define TOFNorm_cxx
#include "TOFNorm.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TOFNorm::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TOFNorm.C
//      root> TOFNorm t
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
    TFile *opf = new TFile("TOFNorm.root", "recreate");
    TTree *opt = new TTree("tree", "tx ntof ce");
    // new data
    Double_t tx, ntof, ce;
    // new tree
    opt->Branch("tx", &tx, "tx/D");
    fChain->AddFriend(opt);

    // itu-itd
    opt->SetAlias("tLimit", "itu < 4095 && itd < 4095 && itu >= 0 && itd >= 0");
    opt->SetAlias("vtLimit", "vitu < 4095 && vitd < 4095 && vitu >= 0 && vitd >= 0");
    fChain->Draw("itd-itu >> tdiff(260, -700, 1900)", "tLimit");


    // dt/dx and fit
    TH1D *tdiff = (TH1D*)gDirectory->Get("tdiff");
    TH1D *dtd = new TH1D("dtd", "dt/dx", 259, -695, 1895);
    Double_t lastBin = tdiff->GetBinContent(1);
    for (int i = 2; i <= tdiff->GetNbinsX(); ++i) {
        Double_t bin = tdiff->GetBinContent(i);
        dtd->Fill(tdiff->GetBinLowEdge(i), bin-lastBin);
        lastBin = bin;
    }
    dtd->Sumw2(0);
    dtd->Draw();
    TF1 *f1 = new TF1("f1", "gaus", -550, -400);
    TF1 *f2 = new TF1("f2", "gaus", 1610, 1720);
    f1->SetParameters(700, -460, 20);
    f2->SetParameters(-600, 1650, 20);
    dtd->Fit(f1, "RN", "");
    dtd->Fit(f2, "R", "");
    f1->Draw("same");


    // calibration
    Double_t txl = f1->GetParameter(1);
    Double_t txr = f2->GetParameter(1);
    // tx = 200.0 / (txr - txl) * ((td-tu)-txl) - 100.0
    Double_t tp0 = 200.0 * txl / (txl-txr) - 100.0;
    Double_t tp1 = 200.0 / (txr-txl);
    // cout << tp0 << "    " << tp1 << endl;
    TString stp0(Form("%lf", tp0));
    TString stp1(Form("%.4lf", tp1));
    // cout << tp0 + tp1 * -700 << "  " << tp0 + tp1 * 1900 << endl;
    fChain->Draw(stp0+"+"+stp1+"*(itd-itu) >> htx(520, -121.25, 121.629)", "tLimit");
    delete f1, f2;

    Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
        fChain->GetEntry(jentry);
        tx = tp0 + tp1 * (itd - itu);
        opt->Fill();
    }

    // opt->Print();

    // test the code before
    fChain->Draw("x:tx-x >> (240, -12, 12, 240, -120, 120)", "tLimit", "colz");



    // fit to fix the tof and get ntof
    // see the ctof
    fChain->Draw("(itu+itd):tx", "tLimit", "colz");

    // profile and fit
    fChain->Draw("(itu+itd)/2.0/40.0:tx >> ctofx", "tLimit && (itu+itd)<2400", "goff");
    TH2D *ctofx = (TH2D*)gDirectory->Get("ctofx");
    TProfile *ctofxpx = ctofx->ProfileX();
    ctofxpx->GetYaxis()->SetRangeUser(27.8, 28.6);
    ctofxpx->Draw();
    f1 = new TF1("f1", "-[0]+[1]*sqrt(x*x+502.5*502.5)/100.0", -110, 110);
    f1->SetParName(0, "TOF_fix");
    f1->SetParName(1, "ntof");
    ctofxpx->Fit(f1, "R");
    Double_t tofFix = f1->GetParameter(0);
    delete f1;


    TBranch *ntofBranch = opt->Branch("ntof", &ntof, "ntof/D");
    TBranch *ceBranch = opt->Branch("ce", &ce, "ce/D");

    for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
        fChain->GetEntry(jentry);
        // ntof
        Double_t ctof = (itu + itd) / 80.0;
        ntof = (ctof + tofFix) * 100.0 / sqrt(502.5*502.5 + tx*tx);
        // ce
        if (ntof > 5) {
            ce = 72.29824 * 72.29824 / ntof / ntof;
        } else {
            ce = 0.0;
        }
        ntofBranch->Fill();
        ceBranch->Fill();
    }

    opt->Print();

    opt->Write();
    opf->Close();
}

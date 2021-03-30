#define AddVeto_cxx
#include "AddVeto.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AddVeto::Loop()
{
//   In a ROOT session, you can do:
//      root> .L AddVeto.C
//      root> AddVeto t
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


    // new file and new tree
    TFile *opf = new TFile("AddVeto.root", "recreate");
    TTree *opt = new TTree("tree", "ADC tree with veto wall");


    // constant for new scintillation
    const Double_t D = 490.0;           // cm, 10cm before the old one
    const Double_t L = 100.0;           // cm, half length
    const Double_t dD = 1.0;            // cm, thickness
    const Double_t TRes = 1.0;          // ns
    const Double_t Lambda = 380.0;      // cm
    const Double_t QRes = 0.1;
    const Double_t Vsc = 7.5;
    // proton
    const Double_t Ec0 = 50.0;          // MeV
    const Double_t EcRes = 50.0;        // MeV
    // ADC
    const Double_t ADCgain = 60.0;
    const Double_t ADCuPed = 120.0;
    const Double_t ADCdPed = 140.0;
    const Double_t ADCnoise = 50.0;
    const Int_t ADCoverflow = 4095;
    // TDC
    const Double_t TriggerDelay = 20.0;
    const Double_t TDCch2ns = 40.0;
    const Int_t TDCoverflow = 4095;
    const Double_t tu_off = 6.4;
    const Double_t td_off = 15.8;


    // variable for new scintillation
    Double_t vx;
    Double_t vtof, vtu, vtd;
    Double_t vqu, vqd;
    Int_t vitu, vitd;
    Int_t viqu, viqd;

    // init new tree
    // old data
    opt->Branch("x", &x, "x/D");
    opt->Branch("e", &e, "e/D");            // energy
    opt->Branch("tof", &tof, "tof/D");      // time of flight
    opt->Branch("pid", &pid, "pid/I");
    opt->Branch("tu", &tu, "tu/D");
    opt->Branch("td", &td, "td/D");
    opt->Branch("qu", &qu, "qu/D");
    opt->Branch("qd", &qd, "qd/D");
    opt->Branch("itu", &itu, "itu/I");
    opt->Branch("itd", &itd, "itd/i");
    opt->Branch("iqu", &iqu, "iqu/i");
    opt->Branch("iqd", &iqd, "iqd/i");
    opt->Branch("diff", &diff, "diff/D");
    // new data
    opt->Branch("vx", &vx, "vx/D");
    opt->Branch("vtof", &vtof, "vtof/D");
    opt->Branch("vtu", &vtu, "vtu/D");
    opt->Branch("vtd", &vtd, "vtd/D");
    opt->Branch("vqu", &vqu, "vqu/D");
    opt->Branch("vqd", &vqd, "vqd/D");
    opt->Branch("vitu", &vitu, "vitu/I");
    opt->Branch("vitd", &vitd, "vitd/I");
    opt->Branch("viqu", &viqu, "viqu/I");
    opt->Branch("viqd", &viqd, "viqd/I");

    // random
    TRandom *gr = new TRandom3(0);


    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        // user code
        vx = x * (D+dD/2.0) / 502.5;                    // probable vx
        Double_t Dr = D + gr->Uniform(-0.5, 0.5) * dD;
        Double_t d = TMath::Sqrt(Dr*Dr + vx*vx);        // cm, flight path

        // tof
        if (pid == 2) {         // proton
            vtof = 72.29824 / TMath::Sqrt(e) * d * 0.01;
        } else {
            vtof = -1;
        }

        // tu, td, qu, qd
        if (pid == 2) {
            vtu = vtof + (L-x)/Vsc + gr->Gaus(0, TRes/2.35) + tu_off - TriggerDelay;
            vtd = vtof + (L+x)/Vsc + gr->Gaus(0, TRes/2.35) + td_off - TriggerDelay;
            vtu *= TDCch2ns;
            vtd *= TDCch2ns;

            // energy deposition in the thin plastic
            Double_t q0 = e * 0.1 * ADCgain;
            // resolution
            q0 = gr->Gaus(q0, q0*QRes/2.35);
            vqu = q0 * TMath::Exp(-(L-x)/Lambda);
            vqd = q0 * TMath::Exp(-(L+x)/Lambda);
            // ADC
            vqu += gr->Gaus(ADCuPed, ADCnoise);
            vqd += gr->Gaus(ADCdPed, ADCnoise);
            vqu = vqu < 0 ? 0 : vqu;
            vqd = vqd < 0 ? 0 : vqd;

        } else {
            vtu = TDCoverflow;
            vtd = TDCoverflow;
            vqu = ADCuPed + gr->Gaus(0, ADCnoise);
            vqd = ADCdPed + gr->Gaus(0, ADCnoise);
        }


        // overflow check
        vtu = vtu > TDCoverflow ? TDCoverflow : vtu;
        vtd = vtd > TDCoverflow ? TDCoverflow : vtd;
        vqu = vqu > ADCoverflow ? ADCoverflow : vqu;
        vqd = vqd > ADCoverflow ? ADCoverflow : vqd;

        // digitization
        vitu = Int_t(vtu);
        vitd = Int_t(vtd);
        viqu = Int_t(vqu);
        viqd = Int_t(vqd);

        opt->Fill();
    }

    opf->Write();
    opf->Close();
}

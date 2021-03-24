//将下列代码保存到readTree.cc
//在ROOT环境下 .x readTree.cc


TCanvas *c1;
void readTree() {
    // 1.打开文件，得到TTree指针
    TFile *ipf = new TFile("tree.root");//打开ROOT文件
    if (ipf->IsZombie()) {
        cout << "Error opening file tree.root" << endl;
        exit(-1);
    }
    ipf->cd();
    TTree *tree = (TTree*)ipf->Get("tree");//得到名字为“tree”的TTree指针

    //2. 声明tree的Branch变量

    Double_t x;
    Double_t e;
    int pid;
    Double_t tof, ctof;
    Double_t tu, td;
    Double_t qu, qd;

    //3. 将变量指向对应Branch的地址
    tree->SetBranchAddress("ctof", &ctof);//将ROOT文件内tree内名为"ctof"的branch的数据的指针指向ctof的变量。
    tree->SetBranchAddress("tof", &tof);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("tu", &tu);
    tree->SetBranchAddress("td", &td);
    tree->SetBranchAddress("qu", &qu);
    tree->SetBranchAddress("qd", &qd);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("e", &e);


    //将新数据写入新的ROOT文件 -对应的代码用 ////标出
    //calibration parameters
    Double_t xp0, xp1;
    //new tree parameters
    Double_t tx, qx, ce, ntof;
    TFile *opf = new TFile("tree2.root", "recreate");
    opf->cd();
    TTree *opt = new TTree("tree", "tree");
    opt->Branch("tu", &tu, "tu/D");
    opt->Branch("td", &td, "td/D");
    opt->Branch("qu", &qu, "qu/D");
    opt->Branch("tx", &tx, "tx/D");
    opt->Branch("x", &x, "x/D");
    opt->Branch("e", &e, "e/D");
    opt->Branch("pid", &pid, "pid/I");
    opt->Branch("ctof", &ctof, "ctof/D");
    opt->Branch("tof", &tof, "tof/D");
    opt->Branch("tx", &tx, "tx/D");
    opt->Branch("qx", &qx, "qx/D");
    opt->Branch("ntof", &ntof, "ntof/D");
    opt->Branch("ce", &ce, "ce/D");


    //4. 逐事件读取tree的branch数据
    // c1 = new TCanvas("c1", "c1");
    // c1->cd();
    Long64_t nentries = tree->GetEntries();//得到tree的事件总数

    // 刻度tx
    tree->Draw("td-tu>>tdiff(140, -20, 50)", "", "goff");
    TH1D *tdiff = (TH1D*)gDirectory->Get("tdiff");
    TH1D *dtd = new TH1D("dtd", "dt/dx", 141, -20.25, 50.25);
    Double_t lastBin = tdiff->GetBinContent(1);
    for (int i = 2; i < tdiff->GetNbinsX(); ++i) {
        Double_t bin = tdiff->GetBinContent(i);
        dtd->Fill(tdiff->GetBinLowEdge(i+1), bin-lastBin);
        lastBin = bin;
    }
    dtd->Sumw2(0);
    dtd->Draw();
    TF1 *f1 = new TF1("f1", "gaus", 39, 45);
    TF1 *f2 = new TF1("f2", "gaus", -13, -9);
    f1->SetParameters(-250, 42, 1.0);
    dtd->Fit(f1, "RN", "");
    dtd->Fit(f2, "R", "same");
    f1->Draw("same");
    Double_t txl = f2->GetParameter(1);
    Double_t txr = f1->GetParameter(1);
    // cout << "txl = " << txl << endl << "txr = " << txr << endl;
    xp0 = 200.0 * txl / (txl-txr) - 100.0;
    xp1 = 200.0 / (txr-txl);
    TString sxp0(Form("%lf", xp0));
    TString sxp1(Form("%lf", xp1));
    tree->Draw(sxp0+"+"+sxp1+"*"+"(td-tu) >> htx(500, -125, 125)");


    for (Long64_t jentry = 0; jentry < nentries; jentry++) {//对事件进行遍历
        tree->GetEntry(jentry);//将第jentry个事件数据填入对应变量(步骤3.中指向的变量)，每次变量值会变成当前事件对应的数据。
        // calculate new parameters
        tx = td - tu;
        tx = xp0 + xp1 * tx;


        //opt->Fill();//fill new parameter to TTree* opt

        // if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
    }


    ipf->Close();
    // opt->Write();
    // opf->Close();
}
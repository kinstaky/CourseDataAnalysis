//将下列代码保存到readTree.cc
//在ROOT环境下 .x readTree.cc

// 物理量定义
const Double_t c = 2.998e8;                     // 光速m/s
const Double_t nc = c * 1e-9;                   // 光速m/ns
const Double_t realTof0 = 5.025 / c * 1e9;      // 光子打在502.5cm外的探测器需要的飞行时间
const Double_t Mn = 939.57;                     // MeV/c^2


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


    //calibration parameters
    Double_t tp0, tp1;
    Double_t txl, txr;
    Double_t qp0, qp1;
    Double_t qxl, qxr;
    //new four parameters
    Double_t tx, qx, ce, ntof;
    TFile *opf = new TFile("tree2.root", "recreate");
    opf->cd();
    TTree *opt = new TTree("tree", "tree");
    // old data
    opt->Branch("tu", &tu, "tu/D");
    opt->Branch("td", &td, "td/D");
    opt->Branch("qu", &qu, "qu/D");
    opt->Branch("qd", &qd, "qd/D");
    opt->Branch("x", &x, "x/D");
    opt->Branch("e", &e, "e/D");
    opt->Branch("pid", &pid, "pid/I");
    opt->Branch("ctof", &ctof, "ctof/D");
    opt->Branch("tof", &tof, "tof/D");
    // new data
    opt->Branch("tx", &tx, "tx/D");
    opt->Branch("qx", &qx, "qx/D");
    opt->Branch("ntof", &ntof, "ntof/D");
    opt->Branch("ce", &ce, "ce/D");


    // c1 = new TCanvas("c1", "c1");
    // c1->cd();
    Long64_t nentries = tree->GetEntries();//得到tree的事件总数

    // 刻度tx
    tree->Draw("td-tu>>tdiff(280, -20, 50)", "", "");
    // c1->Draw();
    TH1D *tdiff = (TH1D*)gDirectory->Get("tdiff");
    TH1D *dtd = new TH1D("dtd", "dt/dx", 279, -19.875, 49.875);
    Double_t lastBin = tdiff->GetBinContent(1);
    for (int i = 2; i < tdiff->GetNbinsX(); ++i) {
        Double_t bin = tdiff->GetBinContent(i);
        dtd->Fill(tdiff->GetBinLowEdge(i), bin-lastBin);
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
    txl = f2->GetParameter(1);
    txr = f1->GetParameter(1);
    // tx = 200.0 / (txr - txl) * ((td-tu)-txl) - 100.0
    tp0 = 200.0 * txl / (txl-txr) - 100.0;
    tp1 = 200.0 / (txr-txl);
    TString stp0(Form("%5.2lf", tp0));
    TString stp1(Form("%4.2lf", tp1));
    tree->Draw(stp0+"+"+stp1+"*"+"(td-tu) >> htx(500, -125, 125)");
    delete f1, f2;

    // 刻度qx
    tree->Draw("log(qu)-log(qd)>>qdiff(320, -0.8, 0.8)", "", "");
    TH1D *qdiff = (TH1D*)gDirectory->Get("qdiff");
    TH1D *dqd = new TH1D("dqd", "dq/dx", 319, -0.7975, 0.7975);
    lastBin = qdiff->GetBinContent(1);
    for (int i = 2; i < qdiff->GetNbinsX(); ++i) {
        Double_t bin = qdiff->GetBinContent(i);
        dqd->Fill(qdiff->GetBinLowEdge(i), bin-lastBin);
        lastBin = bin;
    }
    dqd->Sumw2(0);
    dqd->Draw();
    f1 = new TF1("f1", "gaus", -0.7, -0.3);
    f2 = new TF1("f2", "gaus", 0.35, 0.8);
    f1->SetParameters(300, -0.45, 0.2);
    f2->SetParameters(-300, 0.55, 0.2);
    dqd->Fit(f1, "RN", "");
    dqd->Fit(f2, "R", "same");
    f1->Draw("same");
    qxl = f1->GetParameter(1);
    qxr = f2->GetParameter(1);
    // qx = 200.0 / (qxr - qxl) * ((log(qu)-log(qd))-qxl) - 100.0
    qp0 = 200.0 * qxl / (qxl - qxr) - 100.0;
    qp1 = 200.0 / (qxr-qxl);
    TString sqp0(Form("%5.2lf", qp0));
    TString sqp1(Form("%5.2lf", qp1));
    tree->Draw(sqp0+"+"+sqp1+"*"+"(log(qu)-log(qd)) >> hqx(600, -150, 150)");
    delete f1, f2;

    // 将数据tx和qx暂时填入xTree中
    TTree* xTree = new TTree("xTree", "temporary tree to store tx and qx");
    xTree->Branch("tx", &tx, "tx/D");
    xTree->Branch("qx", &qx, "qx/D");
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {    //对事件进行遍历
        tree->GetEntry(jentry);     //将第jentry个事件数据填入对应变量(步骤3.中指向的变量)，每次变量值会变成当前事件对应的数据。
        // calculate new parameters
        // tx
        tx = td - tu;
        tx = tp0 + tp1 * tx;
        // qx
        qx = log(qu) - log(qd);
        qx = qp0 + qp1 * qx;

        xTree->Fill();
    }
    tree->AddFriend(xTree);

    // 分析一下tx和qx
    tree->Draw("x:tx-x >> h2tx(120, -12, 12, 120, -120, 120)", "", "colz");
    TH2D *h2tx = (TH2D*)(gDirectory->Get("h2tx"));
    TH1D *htxRes = h2tx->ProjectionX();
    htxRes->Fit("gaus");

    tree->Draw("x:qx-x>>h2qx(120, -60, 60, 120, -120, 120)", "", "colz");
    TH2D *h2qx = (TH2D*)(gDirectory->Get("h2qx"));
    TH1D *hqxRes = h2qx->ProjectionX();
    hqxRes->Fit("gaus");



    // 修正TOF
    // 绝对刻度方法
    tree->Draw("ctof >> hgctof(100, 39, 48)", "abs(tx)<5 && ctof<45", "goff");
    TH1D *hgctof = (TH1D*)gDirectory->Get("hgctof");
    f1 = new TF1("f1", "gaus", 42, 45);
    hgctof->Fit(f1, "R");
    Double_t tof0 = f1->GetParameter(1);
    Double_t tofFix = realTof0 - tof0;
    delete f1;

    // 拟合方法
    tree->Draw("ctof:tx >> ctofx", "ctof < 45", "goff");
    TH2D *ctofx = (TH2D*)gDirectory->Get("ctofx");
    TProfile *ctofxpx = ctofx->ProfileX();
    ctofxpx->Sumw2(0);
    ctofxpx->GetYaxis()->SetRangeUser(42.8, 43.5);
    ctofxpx->Draw();
    f1 = new TF1("f1", "-[0]+[1]*sqrt(x*x+502.5*502.5)/100.0", -100, 100);
    f1->SetParName(0, "TOF_fix");
    f1->SetParName(1, "ntof");
    ctofxpx->Fit(f1, "R");
    c1->Draw();
    delete f1;



    for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
        // tx and qx in xTree also getted as friend of tree
        tree->GetEntry(jentry);

        // 计算ntof，归一化到100cm
        ntof = (ctof + tofFix) * 100.0 / sqrt(502.5*502.5 + tx*tx);
        // 计算中子能量ce，考虑了相对论效应
        if (ctof > 45) {
            ce = Mn * (1.0 / sqrt(1.0 - 1.0 / (nc*nc*ntof*ntof)) - 1);
        } else {
            ce = 0.0;
        }
        opt->Fill();
    }

    // 查看新生成的树
    opt->Print();
    // 查看归一化后的光子的飞行时间
    opt->Draw("ntof:tx >> hgctofx(100, -120, 120, 100, 2.8, 3.8)", "", "colz");
    // 查看归一化后的光子和中子的飞行时间
    opt->Draw("ntof:tx >> hctofx(100, -120, 120, 100, 2, 20)", "", "colz");
    // 查看中子能量
    opt->Draw("ce", "ce > 0.0", "");


    ipf->Close();
    opt->Write();
    opf->Close();
}
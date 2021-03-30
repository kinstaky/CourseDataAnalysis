// 将下列代码保存到treeADC.cc
// 在ROOT环境中运行：
// .x treeADC.cc
// 生成 treeADC.root 文件
void treeADC() {
	const Double_t D = 500.0;			// cm, distance between target and the scin.(Center)
	const Double_t L = 100.0;			// cm, half length of the scin.
	const Double_t dD = 5.0;			// cm, thickness of the scin.
	const Double_t TRes = 1.0;			// ns, time resolution(FWHM) of the scintillator.
	const Double_t Lambda = 380.0;		// cm, attenuation lenght of the scin.
	const Double_t QRes = 0.1;			// relative energy resolution(FWHM) of the scin.
	const Double_t Vsc = 7.5;			// ns/cm, speed of light in the scin.
	// neutron & gamma
	const Double_t En0 = 50.0;			// MeV, average neutron energy
	const Double_t EnRes = 50.0;		// MeV, energy spread of neutron(FWHM)
	const Double_t Eg0 = 10.0;			// MeV, gamma energy

	const Double_t Rn = 0.75;			// ratio of neutron
	const Double_t Rg = 0.05;			// ratio of gamma;

	// charged particle hited in Neutron Detector,
	const Double_t Rc = 0.1;			// ratio of proton
	const Double_t Ec0 = 50.0;			// MeV
	const Double_t EcRes = 50.0;		// MeV energy spread of proton(FWHM)

	// ADC
	const Double_t ADCgain = 60.0;		// 1MeV = 60ch
	const Double_t ADCuPed = 140.0;		// baseline of ADC of upper side
	const Double_t ADCdPed = 130.0;		// baseline of ADC of bottom side
	const Double_t ADCnoise = 50.0;		// sigma of noise
	const Int_t ADCoverflow = 4095;

	// time offset
	const Double_t tu_off = 5.5;
	const Double_t td_off = 20.4;

	// TDC
	const Double_t TriggerDelay = 15.0;	// ns, trigger延迟,将感兴趣的时间信号放在TDC量程以内。
	const Double_t TDCch2ns = 40.0;		// 1ns = 40ch
	const Int_t TDCoverflow = 4095;

	TFile *opf = new TFile("treeADC.root", "recreate");		// new file
	TTree *opt = new TTree("tree", "tree structure");		// new tree

	// 定义tree的branch变量
	Double_t x;
	Double_t e;
	Int_t pid;    			// 0/1/2/3: gamma/neutron/proton/No_particle
	Double_t tof, ctof;
	Double_t tu, td;
	Double_t qu, qd;
	Int_t itu,itd;			// TDC
	Int_t iqu,iqd;			// ADC
	Double_t diff;


	opt->Branch("x", &x, "x/D");
	opt->Branch("e", &e, "e/D");			// energy
	opt->Branch("tof", &tof, "tof/D");		// time of flight
	// opt->Branch("ctof",&ctof, "ctof/D");	// TOF from exp. data
	opt->Branch("pid", &pid, "pid/I");
	// raw time and energy
	opt->Branch("tu", &tu, "tu/D");
	opt->Branch("td", &td, "td/D");
	opt->Branch("qu", &qu, "qu/D");
	opt->Branch("qd", &qd, "qd/D");

	// energy in ADC, time in TDC 注意，以下Branch变量声明的类型为 Integer，
	opt->Branch("itu", &itu, "itu/I");
	opt->Branch("itd", &itd, "itd/I");
	opt->Branch("iqu", &iqu, "iqu/I");
	opt->Branch("iqd", &iqd, "iqd/I");

	opt->Branch("diff", &diff, "diff/D");

	TRandom3 *gr = new TRandom3(0);

	for (int i = 0; i < 1000000; ++i) {
		x = gr->Uniform(-L, L);
		Double_t Dr = D + gr->Uniform(-0.5, 0.5) * dD;
		Double_t d = TMath::Sqrt(Dr*Dr + x*x);	// cm, flight path
		// pid
		Double_t ration = gr->Uniform();
		if (ration < Rg) {
			pid = 0;
		} else if (ration < Rg+Rn) {
			pid = 1;
		} else if (ration < Rg+Rn+Rc) {
			pid = 2;
		} else {
			pid = 3;
		}

		// energy & tof
		switch (pid) {
			case 0:				// gamma
				e = Eg0;
				tof = 3.333 * d * 0.01;
				break;
			case 1:				// neutron
				e = gr->Gaus(En0, EnRes/2.35);
				tof = 72.29824 / TMath::Sqrt(e) * d * 0.01;
				break;
			case 2:				// proton
				e = gr->Gaus(Ec0, EcRes/2.35);
				tof = 72.29824 / TMath::Sqrt(e) * d * 0.01;
				break;
			case 3:
				e = -1;
				tof = -1;
				break;
			default:
				cout << "Error: Undefined pid." << endl;
				exit(-1);
		}

		// tu, td, qu, qd
		if (pid == 3) {
			tu = TDCoverflow;
			td = TDCoverflow;
			qu = ADCuPed + gr->Gaus(0, ADCnoise);
			qd = ADCdPed + gr->Gaus(0, ADCnoise);
		} else {
			// time
			tu = tof + (L-x)/Vsc + gr->Gaus(0, TRes/2.35) + tu_off - TriggerDelay;
			td = tof + (L+x)/Vsc + gr->Gaus(0, TRes/2.35) + td_off - TriggerDelay;
			tu *= TDCch2ns;
			td *= TDCch2ns;

			// energy deposition in the plastic
			Double_t q0;
			if (pid != 2) {				// gamma & neutron
				q0 = e * ADCgain * gr->Uniform();
			} else {					// charged particle, full energy
				q0 = e * ADCgain;
			}

			// resolution
			q0 = gr->Gaus(q0, q0*QRes/2.35);

			// light transmissino within the plastic
			qu = q0 * TMath::Exp(-(L-x)/Lambda);
			qd = q0 * TMath::Exp(-(L+x)/Lambda);

			// ADC
			qu += gr->Gaus(ADCuPed, ADCnoise);
			qd += gr->Gaus(ADCdPed, ADCnoise);
			qu = qu < 0 ? 0 : qu;
			qd = qd < 0 ? 0 : qd;
		}

		// overflow check
		tu = tu > TDCoverflow ? TDCoverflow : tu;
		td = td > TDCoverflow ? TDCoverflow : td;
		qu = qu > ADCoverflow ? ADCoverflow : qu;
		qd = qd > ADCoverflow ? ADCoverflow : qd;

		// digitization
		itu = Int_t(tu);
		itd = Int_t(td);
		iqu = Int_t(qu);
		iqd = Int_t(qd);

		// difference in amplitude before and after digitization
		diff = tu - itu;

		opt->Fill();
	}

	opt->Write();
	opf->Close();
}
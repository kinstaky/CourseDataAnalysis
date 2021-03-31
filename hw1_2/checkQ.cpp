void checkQ() {
	gROOT->ProcessLine(".L QProcessor.cpp");
    gROOT->ProcessLine("QProcessor qp");
    gROOT->ProcessLine("qp.Loop()");
}
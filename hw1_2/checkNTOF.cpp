void checkNTOF() {
	gROOT->ProcessLine(".L TOFNorm.cpp");
    gROOT->ProcessLine("TOFNorm t");
    gROOT->ProcessLine("t.Loop()");
}
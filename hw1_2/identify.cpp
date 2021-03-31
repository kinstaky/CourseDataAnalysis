void identify() {
	TFile *ipf = new TFile("AddVeto.root");
	TTree *tree = (TTree*)ipf->Get("tree");
	tree->AddFriend("tree", "TOFNorm.root");
	tree->AddFriend("tree", "QProcessor.root");

	//tree->Draw("q0:ntof >> (480, 4, 20, 450, 0, 4500) ", "qLimit && tLimit && ntof > 5", "colz");

	tree->Draw("q0:ce >> (480, 0, 160, 450, 0, 4500)", "tLimit && qLimit && ce > 0", "colz");

	tree->Draw("q0:ce >> (480, 0, 160, 450, 0, 4500)", "tLimit && qLimit && ce > 0 && !vtLimit", "colz");

	ipf->Close();
}
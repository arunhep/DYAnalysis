{

	const char *histo[2] = {"h_preFSRMass", "h_postFSRMass"};
	const char *titles[2] = {"preFSR Mass (GeV)", "postFSR Mass (GeV)"};
	const char *titles1[2] = {"", ""};
	char name[200];

	TFile *f1 = TFile::Open("DYEE_M10to50_new.root");
	TFile *f2 = TFile::Open("DYEE_M50to120_new.root"); 
	TFile *f3 = TFile::Open("DYEE_M120to200_new.root");
	TFile *f4 = TFile::Open("DYEE_M200to400_new.root");
	TFile *f5 = TFile::Open("DYEE_M400to800_new.root");
	TFile *f6 = TFile::Open("DYEE_M800to1400_new.root");
	TFile *f7 = TFile::Open("DYEE_M1400to2300_new.root"); 
	TFile *f8 = TFile::Open("DYEE_M2300to3500_new.root");

	for(int i=0;i<1;i++){

		//    sprintf(name,histo[i]);

		//**********************************************Z Mass************************
		TCanvas *c1 =  new TCanvas("c1","",200,150,800,550);
		c1->Draw();
		c1->cd();
		//c1->SetLogx();
		//c1->SetLogy();
		c1->SetTickx(1);
		c1->SetTicky(1);

		TH1D *h1 = (TH1D*)f1->Get(histo[i]);
		TH1D *h2 = (TH1D*)f2->Get(histo[i]);
		TH1D *h3 = (TH1D*)f3->Get(histo[i]);
		TH1D *h4 = (TH1D*)f4->Get(histo[i]);
		TH1D *h5 = (TH1D*)f5->Get(histo[i]);
		TH1D *h6 = (TH1D*)f6->Get(histo[i]);
		TH1D *h7 = (TH1D*)f7->Get(histo[i]);
		TH1D *h8  = (TH1D*)f8->Get(histo[i]);

		h1->SetLineColor(kBlack);
		h1->SetFillColor(kBlack);
		h2->SetLineColor(kBlack+1);
		h2->SetFillColor(kBlack+1);
		h3->SetLineColor(kBlack+2);
		h3->SetFillColor(kBlack+2);
		h4->SetLineColor(kBlack+3);
		h4->SetFillColor(kBlack+3);
		h5->SetLineColor(kBlack+4);
		h5->SetFillColor(kBlack+4);
		h6->SetLineColor(kBlack+5);
		h6->SetFillColor(kBlack+5);
		h7->SetLineColor(kBlack+6);
		h7->SetFillColor(kBlack+6);
		h8->SetLineColor(kBlack+7);
		h8->SetFillColor(kBlack+7);

		h1->SetStats(0);
		h2->SetStats(0);
		h3->SetStats(0);
		h4->SetStats(0);
		h5->SetStats(0);
		h6->SetStats(0);
		h7->SetStats(0);
		h8->SetStats(0);

		THStack *hStack1 = new THStack("hs","");

		hStack1->Add(h8);
		hStack1->Add(h7);
		hStack1->Add(h6);
		hStack1->Add(h5);
		hStack1->Add(h4);
		hStack1->Add(h3);
		hStack1->Add(h2);
		hStack1->Add(h1);

		hStack1->Draw("hist");
		hStack1->SetTitle("");
		hStack1->GetXaxis()->SetRangeUser(45.,55.);
		hStack1->GetXaxis()->SetMoreLogLabels();
		hStack1->GetXaxis()->SetNoExponent();
		hStack1->SetTitle(titles1[i]);
		hStack1->GetXaxis()->SetTitle(titles[i]);

		hStack1->GetYaxis()->SetTitle("Number of Events");
		hStack1->SetMaximum(200000);
		//hStack1->GetYaxis()->SetTitleOffset(0.90);
		hStack1->GetYaxis()->SetTitleOffset(1.25);
		hStack1->GetYaxis()->SetTitleSize(0.04);
		hStack1->SetMinimum(0.001);
		hStack1->SetMaximum(5000);

		TLegend *leg1 = new TLegend(0.6905444,0.50,0.8825215,0.85,NULL,"brNDC");
		leg1->AddEntry(h1,"DY M10-50","f");
		leg1->AddEntry(h2,"DY M50-120","f");
		//leg1->AddEntry(h3,"DY M120-200","f");
		//leg1->AddEntry(h4,"DY M200-400","f");
		//leg1->AddEntry(h5,"DY M400-800","f");
		//leg1->AddEntry(h6,"DY M800-1400","f");
		//leg1->AddEntry(h7,"DY M1400-2300","f");
		//leg1->AddEntry(h8,"DY M2300-3500","f");

		leg1->SetFillColor(0);
		leg1->SetLineColor(0);
		leg1->SetTextSize(0.02888889);
		leg1->Draw();

		c1->Draw();
		//sprintf(name, "%s.pdf", histo[i]);
		//c1->SaveAs(name);
		//delete c1;
	}//loop

}

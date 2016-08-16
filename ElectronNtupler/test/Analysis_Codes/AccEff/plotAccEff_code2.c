{
//  TFile *f1 = TFile::Open("DY_10_3000_EE_76X_forAccEff.root");
    TFile *f1 = TFile::Open("DY_10to3000_13June.root");

  TH1F *h_AccPass = (TH1F*)f1->Get("h_mass_AccPass");
  TH1F *h_AccTotal = (TH1F*)f1->Get("h_mass_AccTotal");

  TH1F *h_EffTotal = (TH1F*)f1->Get("h_mass_EffTotal");
  TH1F *h_EffPass = (TH1F*)f1->Get("h_mass_EffPass");

  int n = h_AccPass->GetNbinsX();

  for(int i=1;i<n+1; i++){
   double low  = h_AccPass->GetBinLowEdge(i);
   double high = h_AccPass->GetBinLowEdge(i+1);

   double num = h_EffPass->GetBinContent(i);
   double den = h_EffTotal->GetBinContent(i);

  cout << low << "   " << high << "   " << num/den << endl; 
  }


  TMultiGraph *mg1 = new TMultiGraph();
  TGraphAsymmErrors* gr1 = new TGraphAsymmErrors();
  gr1->BayesDivide(h_AccPass,h_AccTotal,"");
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerStyle(23);
  gr1->SetMarkerSize(1.);

  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors();
  gr2->BayesDivide(h_EffPass,h_EffTotal,"");
  gr2->SetMarkerColor(kBlue);
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerSize(1.);

  mg1->SetTitle("Acceptance;postFSR mass M_{ee} [GeV]; Acceptance");
  TCanvas *c1 = new TCanvas("c1","",200,10,600,600);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogx();
  mg1->Add(gr1);
  mg1->Add(gr2);
  mg1->Draw("APE");
  mg1->GetXaxis()->SetTitleOffset(1.2);
  mg1->GetXaxis()->SetMoreLogLabels();
  mg1->GetXaxis()->SetNoExponent();
  mg1->GetXaxis()->SetRangeUser(10,3000);
  mg1->GetYaxis()->SetTitleOffset(1.4);
  mg1->GetYaxis()->SetRangeUser(0.0,1.015);

  TLegend *leg1 = new TLegend(0.45,0.14,0.65,0.24,NULL,"brNDC");
  leg1->AddEntry(gr1,"Acceptance","lp");
  leg1->AddEntry(gr2,"Event Eff","lp");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.02888889);
  leg1->Draw();
  c1->SaveAs("Acceptance_newMethod_76X_13June.png");

}

{

  TFile *f1 = TFile::Open("SingleElectron_ScaleCorr_new.root"); //Data
  TFile *f2 = TFile::Open("DYEE_M10to3000_new.root"); //DYEE
  TFile *f3 = TFile::Open("DYTT_M10to3000_new.root"); //DYTauTau
  TFile *f4 = TFile::Open("dieleBkg_2_new.root"); //TTbar
  TFile *f5 = TFile::Open("dieleBkg_3_new.root"); //WW
  TFile *f6 = TFile::Open("dieleBkg_4_new.root"); //WZ
  TFile *f7 = TFile::Open("dieleBkg_5_new.root"); //ZZ
  TFile *f8 = TFile::Open("WJets_new.root"); //WJets


  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1","",200,150,800,550);
  c1->Draw();
  c1->cd();
  
 TPad *c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetLogy();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  c1_2->SetGridx();
  c1_2->SetGridy();

  TH1F *h1 = (TH1F*)f1->Get("dieleMass_Peak");
  TH1F *h2 = (TH1F*)f2->Get("dieleMass_Peak");
  TH1F *h3 = (TH1F*)f3->Get("dieleMass_Peak");
  TH1F *h4 = (TH1F*)f4->Get("dieleMass_Peak");
  TH1F *h5 = (TH1F*)f5->Get("dieleMass_Peak");
  TH1F *h6 = (TH1F*)f6->Get("dieleMass_Peak");
  TH1F *h7 = (TH1F*)f7->Get("dieleMass_Peak");
  TH1F *h8  = (TH1F*)f8->Get("dieleMass_Peak");

  TH1F *h_diBoson = (TH1F*)h5->Clone("h_diBoson");
  h_diBoson->Add(h6);
  h_diBoson->Add(h7);

  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8); 
  h1->SetLineColor(1);

  h2->SetLineColor(kOrange-2);
  h2->SetFillColor(kOrange-2);
  h3->SetLineColor(2);  //Purple
  h3->SetFillColor(2);
  h4->SetLineColor(kOrange-8);
  h4->SetFillColor(kOrange-8);
  h_diBoson->SetLineColor(51);
  h_diBoson->SetFillColor(51);
  h8->SetLineColor(kGreen-2);
  h8->SetFillColor(kGreen-2);

  h1->SetStats(0);
  h2->SetStats(0);

  THStack *hStack1 = new THStack("hs","");
  hStack1->Add(h8);
  hStack1->Add(h_diBoson);
  hStack1->Add(h4);
  hStack1->Add(h3);
  hStack1->Add(h2);

  hStack1->Draw("hist");
  h1->Draw("same p");

  hStack1->GetYaxis()->SetTitle("Number of Events");
  hStack1->GetYaxis()->SetTitleOffset(0.90);
  hStack1->GetYaxis()->SetTitleSize(0.04);
  hStack1->SetMinimum(0.09);
  hStack1->SetMaximum(1000000);

  TLegend *leg1 = new TLegend(0.6905444,0.70,0.8825215,0.95,NULL,"brNDC");
  leg1->AddEntry(h1,"data","lp");
  leg1->AddEntry(h2,"dy->ee","f");
  leg1->AddEntry(h3,"dy->#tau#tau","f");
  leg1->AddEntry(h4,"tt(bar)","f");
  leg1->AddEntry(h_diBoson,"diBoson","f");
  leg1->AddEntry(h8,"wjets","f");
  leg1->SetFillColor(0);
  leg1->SetLineColor(0);
  leg1->SetTextSize(0.02888889);
  leg1->Draw();

  c1->cd();

  TH1F *h1_mc = (TH1F *)h2->Clone("h1_mc");
  h1_mc->Add(h3);
  h1_mc->Add(h4);
  h1_mc->Add(h_diBoson);
  h1_mc->Add(h8);

  cout<<"data: "<<h1->Integral()<<endl;
  cout<<"MC: "<<h1_mc->Integral()<<endl;
  cout<<"Ratio:  " << h1->Integral()/h1_mc->Integral() << endl;

  int nbins = h1->GetNbinsX();

  for(int i=1; i<nbins+1; i++){
   double low  = h1->GetBinLowEdge(i);
   double high = h1->GetBinLowEdge(i+1);
//   cout << "Bin  = " << low << "   " << high << endl;
//   cout<<"data: "<<h1->GetBinContent(i)<< endl << "Total MC: "<<h1_mc->GetBinContent(i)<<endl << "Ratio = " << h1->GetBinContent(i)/h1_mc->GetBinContent(i) << endl << "DYEE = " << h2->GetBinContent(i) << endl << "DYTauTau = " << h3->GetBinContent(i) << endl << "ttbar = " << h4->GetBinContent(i) << endl << "diBoson = " << h_diBoson->GetBinContent(i) << endl << "WJets = " << h8->GetBinContent(i) << endl; 
//cout <<"*****************************" << endl;
  }

  TH1F *h1ratio = (TH1F *)h1->Clone("h1ratio");
  h1ratio->Divide(h1_mc);
     
  h1ratio->SetMarkerColor(kBlack);
  h1ratio->SetMarkerStyle(20);
  h1ratio->SetLineColor(kBlack);
  h1ratio->SetMarkerSize(0.8);

  h1ratio->SetTitle("  ");  
  h1ratio->GetXaxis()->SetTitle("Mass (M_{ee}) [GeV]");
  h1ratio->GetXaxis()->SetLabelFont(42);
  h1ratio->GetXaxis()->SetLabelSize(0.15);
  h1ratio->GetXaxis()->SetTitleSize(0.14);
  h1ratio->GetXaxis()->SetTitleOffset(1);

  h1ratio->GetYaxis()->SetTitle("Data/MC");
  h1ratio->GetYaxis()->SetTitleSize(0.15);
  h1ratio->GetYaxis()->SetTitleOffset(0.3);
  h1ratio->GetYaxis()->SetLabelFont(42);
  h1ratio->GetYaxis()->SetLabelSize(0.1);
  h1ratio->GetYaxis()->SetRangeUser(0.4,1.4);
  h1ratio->GetYaxis()->SetNdivisions(5);

  TPad *c1_1 = new TPad("c1_1", "newpad",0.008064516,0.0116071,0.9899194,0.2299107);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetGridx();
  c1_1->SetGridy();
  
  c1_1->Range(-85.9335,-19.83656,785.9335,21.48034);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(1);
  c1_1->SetTopMargin(0.03067478);
  c1_1->SetBottomMargin(0.3047036);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderMode(0);

  h1ratio->Draw();
  TLine l1(15.,1.0,3000.,1.0);
  l1.SetLineWidth(1);
  l1.SetLineColor(kRed);
  l1.SetLineStyle(5);
  l1.Draw("same");
  c1->Draw();
  c1->Update();

  c1->SaveAs("ZMass_Peak.png");

}

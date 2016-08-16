{


  TFile *f1 = TFile::Open("SingleElectron_ScaleCorr_new.root");
  TFile *f2 = TFile::Open("DYEE_M10to3000_new.root");
  TFile *f3 = TFile::Open("DYTT_M10to3000_new.root");
  TFile *f4 = TFile::Open("dieleBkg_2_new.root");
  TFile *f5 = TFile::Open("dieleBkg_3_new.root");
  TFile *f6 = TFile::Open("dieleBkg_4_new.root");
  TFile *f7 = TFile::Open("dieleBkg_5_new.root");
  TFile *f8 = TFile::Open("dieleBkg_1_new.root");

  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1","",323,144,534,523);
  c1->Draw();
  c1->cd();
  
 TPad *c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetLogx();
  c1_2->SetLogy();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
//  c1_2->SetGridx();
  c1_2->SetGridy();

  TH1F *h1_ZMass = (TH1F*)f1->Get("dieleMass_Peak");
  TH1F *h2_ZMass = (TH1F*)f2->Get("dieleMass_Peak");
  TH1F *h3_ZMass = (TH1F*)f3->Get("dieleMass_Peak");
  TH1F *h4_ZMass = (TH1F*)f4->Get("dieleMass_Peak");
  TH1F *h5_ZMass = (TH1F*)f5->Get("dieleMass_Peak");
  TH1F *h6_ZMass = (TH1F*)f6->Get("dieleMass_Peak");
  TH1F *h7_ZMass = (TH1F*)f7->Get("dieleMass_Peak");
  TH1F *h8_ZMass  = (TH1F*)f8->Get("dieleMass_Peak");

  TH1F *h_diBoson_ZMass = (TH1F*)h5_ZMass->Clone("h_diBoson_ZMass");
  h_diBoson_ZMass->Add(h6_ZMass);
  h_diBoson_ZMass->Add(h7_ZMass);

  h1_ZMass->SetMarkerColor(1);
  h1_ZMass->SetMarkerStyle(20);
  h1_ZMass->SetMarkerSize(0.8); 
  h1_ZMass->SetLineColor(1);

  h2_ZMass->SetLineColor(kOrange-2);
  h2_ZMass->SetFillColor(kOrange-2);
  h3_ZMass->SetLineColor(2);  //Purple
  h3_ZMass->SetFillColor(2);
  h4_ZMass->SetLineColor(kOrange-8);
  h4_ZMass->SetFillColor(kOrange-8);
  h_diBoson_ZMass->SetLineColor(51);
  h_diBoson_ZMass->SetFillColor(51);
  h8_ZMass->SetLineColor(kGreen-2);
  h8_ZMass->SetFillColor(kGreen-2);

  h1_ZMass->SetStats(0);
  h2_ZMass->SetStats(0);

  THStack *hStack1 = new THStack("hs","");
  hStack1->Add(h8_ZMass);
  hStack1->Add(h_diBoson_ZMass);
  hStack1->Add(h4_ZMass);
  hStack1->Add(h3_ZMass);
  hStack1->Add(h2_ZMass);

  hStack1->Draw("hist");
  h1_ZMass->Draw("same p");

  hStack1->GetXaxis()->SetMoreLogLabels();
  hStack1->GetXaxis()->SetNoExponent();

  hStack1->GetYaxis()->SetTitle("Number of Events");
  hStack1->GetYaxis()->SetTitleOffset(0.90);
  hStack1->GetYaxis()->SetTitleSize(0.04);
  hStack1->SetMinimum(0.09);
  hStack1->SetMaximum(10000000);

  TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h1_ZMass,"data","lp");
  leg1->AddEntry(h2_ZMass,"dy->ee","f");
  leg1->AddEntry(h3_ZMass,"dy->#tau#tau","f");
  leg1->AddEntry(h4_ZMass,"tt(bar)","f");
  leg1->AddEntry(h_diBoson_ZMass,"diBoson","f");
  leg1->AddEntry(h8_ZMass,"wjets","f");
  leg1->SetFillColor(0);
  leg1->SetLineColor(0);
  leg1->SetTextSize(0.02888889);
  leg1->Draw();

  c1->cd();

  TH1F *h1_mc = (TH1F *)h2_ZMass->Clone("h1_mc");
  h1_mc->Add(h3_ZMass);
  h1_mc->Add(h4_ZMass);
  h1_mc->Add(h_diBoson_ZMass);
  h1_mc->Add(h8_ZMass);

  //cout<<"data: "<<h1_ZMass->Integral()<<endl;
  //cout<<"MC: "<<h1_mc->Integral()<<endl;

  /*int nbins = h1_ZMass->GetNbinsX();
  cout<<"nbins: "<<nbins<<endl;
  cout<<"No Masking"<<endl;

  for(int i=1; i<nbins+1; i++){

   double low  = h1_ZMass->GetBinLowEdge(i);
   double high = h1_ZMass->GetBinLowEdge(i+1);

   cout<<low<<"   "<<high<<"   "<<"data: "<<h1_ZMass->GetBinContent(i)<<"      MC: "<<h1_mc->GetBinContent(i)<<endl;
  }*/

  TH1F *h1ratio = (TH1F *)h1_ZMass->Clone("h1ratio");
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
  h1ratio->GetXaxis()->SetRangeUser(60.,120.);
//  h1ratio->GetXaxis()->SetMoreLogLabels();
//  h1ratio->GetXaxis()->SetNoExponent();

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
//  c1_1->SetLogx();
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
  TLine l1(60.,1.0,120.,1.0);
  l1.SetLineWidth(1);
  l1.SetLineColor(kRed);
  l1.SetLineStyle(5);
  l1.Draw("same");
  c1->Draw();
  c1->Update();

  c1->SaveAs("ZMass_Peak_PUWeighted_76X_ScaleCorr.png");

}

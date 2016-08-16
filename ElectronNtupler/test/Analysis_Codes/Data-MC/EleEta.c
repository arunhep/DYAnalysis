{


  TFile *f1 = TFile::Open("SingleElectron_ScaleCorr_new.root");
  TFile *f2 = TFile::Open("DYEE_M10to3000_All_new.root");
  TFile *f3 = TFile::Open("DYTauTau_M10to3000_All_new.root");
  TFile *f4 = TFile::Open("BKG_2_dzMask_new.root");
  TFile *f5 = TFile::Open("BKG_3_dzMask_new.root");
  TFile *f6 = TFile::Open("BKG_4_dzMask_new.root");
  TFile *f7 = TFile::Open("BKG_5_dzMask_new.root");
  TFile *f8 = TFile::Open("BKG_1_dzMask_new.root");

  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1","",323,144,534,523);
  c1->Draw();
  c1->cd();
  
 TPad *c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
 // c1_2->SetLogx();
  c1_2->SetLogy();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  c1_2->SetGridx();
  c1_2->SetGridy();

  TH1F *h1_ZMass = (TH1F*)f1->Get("eleEta");
  TH1F *h2_ZMass = (TH1F*)f2->Get("eleEta");
  TH1F *h3_ZMass = (TH1F*)f3->Get("eleEta");
  TH1F *h4_ZMass = (TH1F*)f4->Get("eleEta");
  TH1F *h5_ZMass = (TH1F*)f5->Get("eleEta");
  TH1F *h6_ZMass = (TH1F*)f6->Get("eleEta");
  TH1F *h7_ZMass = (TH1F*)f7->Get("eleEta");
  TH1F *h8_ZMass  = (TH1F*)f8->Get("eleEta");

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

//cout << "data = " << h1_ZMass->Integral() << endl;


  hStack1->GetXaxis()->SetMoreLogLabels();
  hStack1->GetXaxis()->SetNoExponent();

  hStack1->GetYaxis()->SetTitle("Number of Events");
  hStack1->GetYaxis()->SetTitleOffset(0.90);
  hStack1->GetYaxis()->SetTitleSize(0.04);
  hStack1->SetMinimum(0.09);
  hStack1->SetMaximum(10000000);

  TLegend *leg1 = new TLegend(0.6905444,0.7559494,0.8825215,0.9517722,NULL,"brNDC");
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

  cout<<"data: "<<h1_ZMass->Integral()<<endl;
  cout<<"MC: "<<h1_mc->Integral()<<endl;

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
  h1ratio->GetXaxis()->SetTitle("Electron #eta");
  h1ratio->GetXaxis()->SetLabelFont(42);
  h1ratio->GetXaxis()->SetLabelSize(0.15);
  h1ratio->GetXaxis()->SetTitleSize(0.14);
  h1ratio->GetXaxis()->SetTitleOffset(1);
//  h1ratio->GetXaxis()->SetRangeUser(15.,3000.);
//  h1ratio->GetXaxis()->SetMoreLogLabels();
//  h1ratio->GetXaxis()->SetNoExponent();

  h1ratio->GetYaxis()->SetTitle("Data/MC");
  h1ratio->GetYaxis()->SetTitleSize(0.15);
  h1ratio->GetYaxis()->SetTitleOffset(0.3);
  h1ratio->GetYaxis()->SetLabelFont(42);
  h1ratio->GetYaxis()->SetLabelSize(0.1);
  h1ratio->GetYaxis()->SetRangeUser(0.6,1.6);
  h1ratio->GetYaxis()->SetNdivisions(8);

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
  TLine l1(15.,1.0,3000.,1.0);
  l1.SetLineWidth(1);
  l1.SetLineColor(kRed);
  l1.SetLineStyle(5);
  l1.Draw("same");
  c1->Draw();
  c1->Update();

  c1->SaveAs("EleEta_mcSmeared_PUWeighted.png");
/*
  // *********************************************Z Mass Peak************************

  TCanvas *c2 =  new TCanvas("c2","",323,144,534,523);
  c2->Draw();
  c2->cd();
  
 TPad *c2_2 = new TPad("c2_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c2_2->SetBottomMargin(0.1);
  c2_2->Draw();
  c2_2->cd();
  //c2_2->SetLogx();
  c2_2->SetLogy();
  c2_2->SetTickx(1);
  c2_2->SetTicky(1);
  c2_2->SetGridx();
  c2_2->SetGridy();

  TH1F *h1_ZMassPeak = (TH1F*)f1->Get("eleEta_Peak");
  TH1F *h2_ZMassPeak = (TH1F*)f2->Get("eleEta_Peak");
  TH1F *h3_ZMassPeak = (TH1F*)f3->Get("eleEta_Peak");
  TH1F *h4_ZMassPeak = (TH1F*)f4->Get("eleEta_Peak");
  TH1F *h5_ZMassPeak = (TH1F*)f5->Get("eleEta_Peak");
  TH1F *h6_ZMassPeak = (TH1F*)f6->Get("eleEta_Peak");
  TH1F *h7_ZMassPeak = (TH1F*)f7->Get("eleEta_Peak");
  TH1F *h8_ZMassPeak = (TH1F*)f8->Get("eleEta_Peak");

  TH1F *h_diBoson = (TH1F*)h5_ZMassPeak->Clone("h_diBoson_ZMassPeak");
  h_diBoson_ZMassPeak->Add(h6_ZMassPeak);
  h_diBoson_ZMassPeak->Add(h7_ZMassPeak);

  h1_ZMassPeak->SetMarkerColor(1);
  h1_ZMassPeak->SetMarkerStyle(20);
  h1_ZMassPeak->SetMarkerSize(0.8); 
  h1_ZMassPeak->SetLineColor(1);

  h2_ZMassPeak->SetLineColor(kOrange-2);
  h2_ZMassPeak->SetFillColor(kOrange-2);
  h3_ZMassPeak->SetLineColor(2);  //Purple
  h3_ZMassPeak->SetFillColor(2);
  h4_ZMassPeak->SetLineColor(kOrange-8);
  h4_ZMassPeak->SetFillColor(kOrange-8);
  h_diBoson_ZMassPeak->SetLineColor(51);
  h_diBoson_ZMassPeak->SetFillColor(51);
  h8_ZMassPeak->SetLineColor(kGreen-2);
  h8_ZMassPeak->SetFillColor(kGreen-2);

  h1_ZMassPeak->SetStats(0);
  h2_ZMassPeak->SetStats(0);

  THStack *hStack2 = new THStack("hs","");
  hStack2->Add(h8_ZMassPeak);
  hStack2->Add(h_diBoson_ZMassPeak);
  hStack2->Add(h4_ZMassPeak);
  hStack2->Add(h3_ZMassPeak);
  hStack2->Add(h2_ZMassPeak);

  hStack2->Draw("hist");
  h1_ZMassPeak->Draw("same p");

  //hStack2->GetXaxis()->SetMoreLogLabels();
  //hStack2->GetXaxis()->SetNoExponent();

  hStack2->GetYaxis()->SetTitle("Number of Events");
  hStack2->GetYaxis()->SetTitleOffset(0.90);
  hStack2->GetYaxis()->SetTitleSize(0.05);
  hStack2->SetMinimum(0.9);
  hStack2->SetMaximum(10000000);

  TLegend *leg2 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg2->AddEntry(h1_ZMassPeak,"data","lp");
  leg2->AddEntry(h2_ZMassPeak,"dy->ee","f");
  leg2->AddEntry(h3_ZMassPeak,"dy->#tau#tau","f");
  leg2->AddEntry(h4_ZMassPeak,"tt(bar)","f");
  leg2->AddEntry(h_diBoson_ZMassPeak,"diBoson","f");
  leg2->AddEntry(h8_ZMassPeak,"wjets","f");
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.02888889);
  leg2->Draw();

  c2->cd();

  TH1F *h2_mc = (TH1F *)h2_ZMassPeak->Clone("h1_mc");
  h2_mc->Add(h3_ZMassPeak);
  h2_mc->Add(h4_ZMassPeak);
  h2_mc->Add(h_diBoson_ZMassPeak);
  h2_mc->Add(h8_ZMassPeak);

  TH1F *h2ratio = (TH1F *)h1_ZMassPeak->Clone("h2ratio");
  h2ratio->Divide(h2_mc);
     
  h2ratio->SetMarkerColor(kBlack);
  h2ratio->SetMarkerStyle(20);
  h2ratio->SetLineColor(kBlack);
  h2ratio->SetMarkerSize(0.8);

  h2ratio->SetTitle("  ");  
  h2ratio->GetXaxis()->SetTitle("Mass (M_{ee})[GeV]");
  h2ratio->GetXaxis()->SetLabelFont(42);
  h2ratio->GetXaxis()->SetLabelSize(0.15);
  h2ratio->GetXaxis()->SetTitleSize(0.14);
  h2ratio->GetXaxis()->SetTitleOffset(1);
  h2ratio->GetXaxis()->SetRangeUser(60.,120.);
  h2ratio->GetXaxis()->SetMoreLogLabels();
  h2ratio->GetXaxis()->SetNoExponent();

  h2ratio->GetYaxis()->SetTitle("Data/MC");
  h2ratio->GetYaxis()->SetTitleSize(0.15);
  h2ratio->GetYaxis()->SetTitleOffset(0.3);
  h2ratio->GetYaxis()->SetLabelFont(42);
  h2ratio->GetYaxis()->SetLabelSize(0.1);
  h2ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  h2ratio->GetYaxis()->SetNdivisions(5);

  TPad *c2_1 = new TPad("c2_1", "newpad",0.008064516,0.0116071,0.9899194,0.2299107);
  c2_1->Draw();
  c2_1->cd();
  //c2_1->SetLogx();
  c2_1->SetGridx();
  c2_1->SetGridy();
  
  c2_1->Range(-85.9335,-19.83656,785.9335,21.48034);
  c2_1->SetFillColor(0);
  c2_1->SetBorderMode(0);
  c2_1->SetBorderSize(1);
  c2_1->SetTopMargin(0.03067478);
  c2_1->SetBottomMargin(0.3047036);
  c2_1->SetFrameBorderMode(0);
  c2_1->SetFrameBorderMode(0);

  h2ratio->Draw();
  TLine l2(60.,1.0,120.,1.0);
  l2.SetLineWidth(1);
  l2.SetLineColor(kRed);
  l2.SetLineStyle(5);
  l2.Draw("same");
  c2->Draw();
  c2->Update();

  c2->SaveAs("ZMassPeak_mcSmeared.png");
  */
}

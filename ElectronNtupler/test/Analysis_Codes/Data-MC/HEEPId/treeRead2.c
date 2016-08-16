{
 const char *file1[1] = {"SingleElectron_UnScaleCorr"};

  char name[200];
  double Ele1PT, Ele2PT, Ele1Eta, Ele2Eta, Ele1Phi, Ele2Phi, ZMass, ZMassCorr, ZPT, ZEta, ZRapidity, ZPhi,primVtx;
  char BB, BE, EE;

  Double_t xbin[46] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1200, 1500, 2000, 3000};
  //Double_t xbin[6] = {830, 1000, 1200, 1500, 2000, 3000};

  TH1F *ele1PT[1]; TH1F *ele2PT[1]; TH1F *elePT[1]; TH1F *ele1Eta[1]; TH1F *ele2Eta[1];
  TH1F *eleEta[1]; TH1F *ele1Phi[1]; TH1F *ele2Phi[1]; TH1F *elePhi[1];
  TH1F *dieleMass[1]; TH1F *dieleMass_BB[1]; TH1F *dieleMass_BE[1]; TH1F *dieleMass_EE[1];
  TH1F *dieleMass_Peak[1];  TH1F *dieleMass_Peak_BB[1]; TH1F *dieleMass_Peak_BE[1]; TH1F *dieleMass_Peak_EE[1];
  TH1F *dielePT[1]; TH1F *dieleEta[1]; TH1F *dieleRap[1]; TH1F *dielePhi[1]; 
  TH1F *nPV[1]; TH1F *nPV_wt[1];

  TFile *file2[1];

  for(int i=0;i<1;i++){

    sprintf(name, "%s.root", file1[i]);
    TFile *f = TFile::Open(name);
    
    sprintf(name, "/afs/cern.ch/work/a/arun/DYAnalysis2015_76X_v1/CMSSW_7_6_3/src/EgammaWork/ElectronNtupler/test/Data-MC/HEEPId/Plots_HEEPId/%s_new.root", file1[i]);
    file2[i] = new TFile(name,"RECREATE");
 
    TTree *t;
    t = (TTree*)f->Get("tree");
    Long64_t entries = t->GetEntries();
    cout<<"entry: "<<entries<<endl;

    t->SetBranchAddress("primVtx",&primVtx);
    t->SetBranchAddress("Ele1PT",&Ele1PT);
    t->SetBranchAddress("Ele2PT",&Ele2PT);
    t->SetBranchAddress("Ele1Eta",&Ele1Eta);
    t->SetBranchAddress("Ele2Eta",&Ele2Eta);
    t->SetBranchAddress("Ele1Phi",&Ele1Phi);
    t->SetBranchAddress("Ele2Phi",&Ele2Phi);
    t->SetBranchAddress("ZMass",&ZMass);
    t->SetBranchAddress("ZPT",&ZPT);
    t->SetBranchAddress("ZEta",&ZEta);
    t->SetBranchAddress("ZRapidity",&ZRapidity);
    t->SetBranchAddress("ZPhi",&ZPhi);
    t->SetBranchAddress("BB",&BB);
    t->SetBranchAddress("BE",&BE);
    t->SetBranchAddress("EE",&EE);

    sprintf(name, "%s",file1[i]);
    nPV[i] = new TH1F("nPV","nPV",40,0.,40.);
    nPV_wt[i] = new TH1F("nPV_wt","nPV Weighted",40,0.,40.);
    ele1PT[i]  = new TH1F("ele1PT", "ele1PT" , 100, 0, 1000);
    ele2PT[i]  = new TH1F("ele2PT", "ele2PT" , 100, 0, 1000);
    elePT[i]   = new TH1F("elePT", "elePT" , 100, 0, 1000);
    ele1Eta[i] = new TH1F("ele1Eta", "ele1Eta" , 50, -2.5, 2.5);
    ele2Eta[i] = new TH1F("ele2Eta", "ele2Eta" , 50, -2.5, 2.5);
    eleEta[i]  = new TH1F("eleEta", "eleEta" , 50, -2.5, 2.5);
    ele1Phi[i] = new TH1F("ele1Phi", "ele1Phi" , 50, -3.5, 3.5);
    ele2Phi[i] = new TH1F("ele2Phi", "ele2Phi" , 50, -3.5, 3.5);
    elePhi[i]  = new TH1F("elePhi", "elePhi" , 50, -3.5, 3.5);
    dieleMass[i] = new TH1F("dieleMass", "dieleMass", 45.0, xbin);
    dieleMass_BB[i] = new TH1F("dieleMass_BB", "dieleMass_BB", 45, xbin);
    dieleMass_BE[i] = new TH1F("dieleMass_BE", "dieleMass_BE", 45, xbin);
    dieleMass_EE[i] = new TH1F("dieleMass_EE", "dieleMass_EE", 45, xbin);
    dieleMass_Peak[i] = new TH1F("dieleMass_Peak", "dieleMass_Peak", 40, 80, 100);
    dieleMass_Peak_BB[i] = new TH1F("dieleMass_Peak_BB", "dieleMass_Peak_BB", 40, 80, 100);
    dieleMass_Peak_BE[i] = new TH1F("dieleMass_Peak_BE", "dieleMass_Peak_BE", 40, 80, 100);
    dieleMass_Peak_EE[i] = new TH1F("dieleMass_Peak_EE", "dieleMass_Peak_EE", 40, 80, 100);
    dielePT[i]   = new TH1F("dielePT", "dielePT" , 100, 0, 1000);
    dieleEta[i]  = new TH1F("dieleEta", "dieleEta" , 60, -8, 8);
    dieleRap[i]  = new TH1F("dieleRap", "dieleRap" , 60, -3, 3);
    dielePhi[i]  = new TH1F("dielePhi", "dielePhi" , 50, -3.5, 3.5);


    ele1PT[i]->Sumw2(); ele1Eta[i]->Sumw2(); ele1Phi[i]->Sumw2();
    ele2PT[i]->Sumw2(); ele2Eta[i]->Sumw2(); ele2Phi[i]->Sumw2();
    elePT[i]->Sumw2(); eleEta[i]->Sumw2(); elePhi[i]->Sumw2();
    dieleMass[i]->Sumw2(); dieleMass_BB[i]->Sumw2(); dieleMass_BE[i]->Sumw2(); dieleMass_EE[i]->Sumw2();
    dieleMass_Peak[i]->Sumw2(); dieleMass_Peak_BB[i]->Sumw2(); dieleMass_Peak_BE[i]->Sumw2(); dieleMass_Peak_EE[i]->Sumw2();
    dielePT[i]->Sumw2(); dieleEta[i]->Sumw2(); dieleRap[i]->Sumw2(); dielePhi[i]->Sumw2();
    nPV[i]->Sumw2();nPV_wt[i]->Sumw2();

    for(Long64_t j=0;j<entries;j++){
      t->GetEntry(j);

      ele1PT[i]->Fill(Ele1PT);
      ele2PT[i]->Fill(Ele2PT);
      elePT[i]->Fill(Ele1PT);
      elePT[i]->Fill(Ele2PT);

      ele1Eta[i]->Fill(Ele1Eta);
      ele2Eta[i]->Fill(Ele2Eta);
      eleEta[i]->Fill(Ele1Eta);
      eleEta[i]->Fill(Ele2Eta);

      ele1Phi[i]->Fill(Ele1Phi);
      ele2Phi[i]->Fill(Ele2Phi);
      elePhi[i]->Fill(Ele1Phi);
      elePhi[i]->Fill(Ele2Phi);

      nPV[i]->Fill(primVtx);
      nPV_wt[i]->Fill(primVtx);
      dieleMass[i]->Fill(ZMass);
      if(BB) dieleMass_BB[i]->Fill(ZMass);
      if(BE) dieleMass_BE[i]->Fill(ZMass);
      if(EE) dieleMass_EE[i]->Fill(ZMass);
      dieleMass_Peak[i]->Fill(ZMass);
      if(BB) dieleMass_Peak_BB[i]->Fill(ZMass);
      if(BE) dieleMass_Peak_BE[i]->Fill(ZMass);
      if(EE) dieleMass_Peak_EE[i]->Fill(ZMass);

      dielePT[i]->Fill(ZPT);
      dieleEta[i]->Fill(ZEta);
      dieleRap[i]->Fill(ZRapidity);
      dielePhi[i]->Fill(ZPhi);
    }

   file2[i]->Write();
   file2[i]->Close();
  }
}


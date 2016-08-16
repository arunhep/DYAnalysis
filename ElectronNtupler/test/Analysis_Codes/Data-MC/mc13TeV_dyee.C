#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t dEta = deltaEta(eta1, eta2);
  Double_t dPhi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(dEta*dEta + dPhi*dPhi);
  return dr;
}

void mc13TeV_dyee() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;

int mass[11] = {10,50,200,400,500,700,800,1000,1500,2000,3000};
//double xsec[10] = {18609.9,6015,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
double xsec[10] = {18609.9,6015,7.68,0.423,0.24,0.036,0.03,0.0159,0.00201,0.00054};


double sumofWts[10] = {2271453516138.643555,451509656182.288330,20960115.865038,368792.490468,209968.687946,32878.675926,27995.408368,14770.905713,1989.393856,475.385433};

//double sumofWts[10] = {2271453516138.643555,451509656182.288330,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

workdir = "/tmp/arun/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";

  InputFiles_signal_DY.clear();

  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_200to400.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_400to500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_500to700.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_700to800.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_800to1000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1000to1500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1500to2000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));

  int nsample = InputFiles_signal_DY.size();
  TFile* file[10]; 

  for(unsigned int jentry = 0; jentry < InputFiles_signal_DY.size(); ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("dataPUDist.root");
    TFile *f2 = TFile::Open("PileUp_MC.root");

    // data histogram 
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Rap;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt; 
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<int>     *passVetoId;
    vector<int>     *passLooseId;
    vector<int>     *passMediumId;
    vector<int>     *passTightId;
    vector<int>     *passHEEPId;
    vector<int>     *isPassMedium_NoPt;
    vector<int>     *isPassMedium_NoScEta;
    vector<int>     *isPassMedium_NoDEta;
    vector<int>     *isPassMedium_NoDPhi;
    vector<int>     *isPassMedium_NoSigmaEtaEta;
    vector<int>     *isPassMedium_NoHOverE;
    vector<int>     *isPassMedium_NoDxy;
    vector<int>     *isPassMedium_NoDz;
    vector<int>     *isPassMedium_NoEInvP;
    vector<int>     *isPassMedium_NoPFIso;
    vector<int>     *isPassMedium_NoConVeto;
    vector<int>     *isPassMedium_NoMissHits;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *rapElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    vector<double>  *pt_Ele23;
    vector<double>  *eta_Ele23;
    vector<double>  *phi_Ele23;
    Int_t           nPV;
    Int_t           nPUTrue;
    vector<float>   *dz;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Rap = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0; 
    genPreFSR_Eta = 0; 
    genPreFSR_Rap = 0; 
    genPreFSR_Phi = 0; 
    genPreFSR_En = 0;
    ptElec = 0;
    etaElec = 0;
    rapElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passVetoId = 0;
    passLooseId = 0;
    passMediumId = 0;
    passTightId = 0;
    passHEEPId = 0;
    isPassMedium_NoPt = 0;
    isPassMedium_NoScEta = 0;
    isPassMedium_NoDEta = 0;
    isPassMedium_NoDPhi = 0;
    isPassMedium_NoSigmaEtaEta = 0;
    isPassMedium_NoHOverE = 0;
    isPassMedium_NoDxy = 0;
    isPassMedium_NoDz = 0;
    isPassMedium_NoEInvP = 0;
    isPassMedium_NoPFIso = 0;
    isPassMedium_NoConVeto = 0;
    isPassMedium_NoMissHits = 0;
    dz = 0;
    pt_Ele23 = 0;
    eta_Ele23 = 0;
    phi_Ele23 = 0;

    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt); 
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta); 
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap); 
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi); 
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("rapElec", &rapElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passVetoId", &passVetoId);
    T1->SetBranchAddress("passLooseId", &passLooseId);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("passTightId", &passTightId);
    T1->SetBranchAddress("passHEEPId", &passHEEPId);
    T1->SetBranchAddress("isPassMedium_NoPt", &isPassMedium_NoPt);
    T1->SetBranchAddress("isPassMedium_NoScEta", &isPassMedium_NoScEta);
    T1->SetBranchAddress("isPassMedium_NoDEta", &isPassMedium_NoDEta);
    T1->SetBranchAddress("isPassMedium_NoDPhi", &isPassMedium_NoDPhi);
    T1->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta);
    T1->SetBranchAddress("isPassMedium_NoHOverE", &isPassMedium_NoHOverE);
    T1->SetBranchAddress("isPassMedium_NoDxy", &isPassMedium_NoDxy);
    T1->SetBranchAddress("isPassMedium_NoDz", &isPassMedium_NoDz);
    T1->SetBranchAddress("isPassMedium_NoEInvP", &isPassMedium_NoEInvP);
    T1->SetBranchAddress("isPassMedium_NoPFIso", &isPassMedium_NoPFIso);
    T1->SetBranchAddress("isPassMedium_NoConVeto", &isPassMedium_NoConVeto);
    T1->SetBranchAddress("isPassMedium_NoMissHits", &isPassMedium_NoMissHits);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("pt_Ele23", &pt_Ele23);
    T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
    T1->SetBranchAddress("phi_Ele23", &phi_Ele23);
    T1->SetBranchAddress("nPV", &nPV);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);
    T1->SetBranchAddress("dz", &dz);

    file[jentry] = new TFile(Form("MediumId_Calibrated/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double massGen;
    int mediumId, dzMaskId, passId, heepId;
    bool trigMatch_lead, trigMatch_slead;
    double dR, dR1, dR2;
    double sum_weights;
    TLorentzVector ele1,ele2,dielectron;
    TLorentzVector gen1,gen2,diGen;

    // Branch variable declaration
    double primVtx;
    double Ele1PT, Ele2PT, Ele1Eta, Ele2Eta, Ele1Phi, Ele2Phi;
    double ZMass, ZPT, ZEta, ZRapidity, ZPhi, DZ1, DZ2;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE;

    // Branch declaration
    tree->Branch("primVtx", &primVtx, "primVtx/D");
    tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
    tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
    tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
    tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
    tree->Branch("Ele1Phi", &Ele1Phi, "Ele1Phi/D");
    tree->Branch("Ele2Phi", &Ele2Phi, "Ele2Phi/D");
    tree->Branch("ZMass", &ZMass, "ZMass/D");
    tree->Branch("ZPT", &ZPT, "ZPT/D");
    tree->Branch("ZEta", &ZEta, "ZEta/D");
    tree->Branch("ZRapidity", &ZRapidity, "ZRapidity/D");
    tree->Branch("ZPhi", &ZPhi, "ZPhi/D");
    tree->Branch("DZ1", &DZ1, "DZ1/D");
    tree->Branch("DZ2", &DZ2, "DZ2/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    tree->Branch("BB", &BB, "BB/B");
    tree->Branch("BE", &BE, "BE/B");
    tree->Branch("EE", &EE, "EE/B");

    vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; 
    vector <double> newscEta; vector <double> neweleDZ;

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    sum_weights = 0.;

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      int index[ptElec->size()];
      float pt[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++)
      {
	pt[el]=ptElec->at(el);
      }

      int size = sizeof(pt)/sizeof(pt[0]);
      TMath::Sort(size,pt,index,true);

      trigMatch_lead = false; trigMatch_slead = false;
      dR = 0.; dR1 = 0.; dR2 = 0.;
      BB=false; BE=false; EE=false;
      mediumId = 0; dzMaskId = 0; passId = 0; heepId = 0;
      newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); newscEta.clear(); neweleDZ.clear();

      // PU Weight
      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);
      //cout << "True PU = " << nPUTrue << "   " << "bin =  " << bin << "  " << "PUWeight =  " << PUWeight << endl;

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();

	if(!tauFlag) (sum_weights = sum_weights + theWeight);

      }

      if(jentry==1){if(massGen > 200.) continue;}        // Gen Mass cut ----- for 50 to inf sample
  //    if(jentry==5){if(massGen < 700. || massGen > 800.) continue;}        // Gen Mass cut ----- for 50 to inf sample
  //    if(jentry==6){if(massGen < 800. || massGen > 1000.) continue;}        // Gen Mass cut ----- for 50 to inf sample
  //    if(jentry==7){if(massGen < 1000. || massGen > 1500.) continue;}        // Gen Mass cut ----- for 50 to inf sample
  //    if(jentry==8){if(massGen < 1500. || massGen > 2000.) continue;}        // Gen Mass cut ----- for 50 to inf sample
  //    if(jentry==9){if(massGen < 2000. || massGen > 3000.) continue;}        // Gen Mass cut ----- for 50 to inf sample


      if(tauFlag) continue;

      if(!Ele23_WPLoose) continue;

      primVtx = nPV;

      if(ptElec->size()>=2.){

	for(int j=0;j<ptElec->size();j++){

//	  heepId = passHEEPId->at(index[j]);
	  mediumId = passMediumId->at(index[j]);
	  //dzMaskId = isPassMedium_NoDz->at(index[j]);

	  //if (ptElec->at(j) < 500.) { passId = mediumId; }
	  //else { passId = dzMaskId; }

	  if(mediumId) {

	    if(fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566)){

	      newelePt.push_back(ptElec->at(index[j]));
	      neweleEta.push_back(etaElec->at(index[j]));
	      neweleEnr.push_back(energyElec->at(index[j]));
	      newelePhi.push_back(phiElec->at(index[j]));
	      newscEta.push_back(etaSC->at(index[j]));
	      neweleDZ.push_back(dz->at(index[j]));


	    } // eta
	  } // ID
	} // pt size
      } // size >=2

      if(newelePt.size()==2){

	for(unsigned int j = 0; j < pt_Ele23->size(); j++){

	  double dR1_comp = 1000.;
	  double dR2_comp = 1000.;

	  dR1 = deltaR(neweleEta.at(0), newelePhi.at(0), eta_Ele23->at(j), phi_Ele23->at(j));
	  dR2 = deltaR(neweleEta.at(1), newelePhi.at(1), eta_Ele23->at(j), phi_Ele23->at(j));

	  if(dR1 < 0.1){
	    if (dR1 < dR1_comp)
	    {
	      dR1_comp = dR1;
	      trigMatch_lead = true;
	    }
	  } // dR1

	  if(dR2 < 0.1){
	    if (dR2 < dR2_comp)
	    {
	      dR2_comp = dR2;
	      trigMatch_slead = true;
	    }
	  } // dR2
	} // pt_Ele23

	if(trigMatch_lead || trigMatch_slead){

	  if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	    Ele1PT  = newelePt.at(0);
	    Ele2PT  = newelePt.at(1);
	    Ele1Eta = neweleEta.at(0);
	    Ele2Eta = neweleEta.at(1);
	    Ele1Phi = newelePhi.at(0);
	    Ele2Phi = newelePhi.at(1);
	    DZ1 = neweleDZ.at(0);
	    DZ2 = neweleDZ.at(1);

	    ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));
	    dielectron=ele1+ele2;

	    ZMass = dielectron.M();
	    ZPT = dielectron.Pt();
	    ZEta = dielectron.Eta();
	    ZRapidity = dielectron.Rapidity();
	    ZPhi = dielectron.Phi();

	    lumiWeight = lumi_Weight;
	    genWeight  = theWeight; 

	    if(fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) < 1.4442) BB = true;
	    if((fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) > 1.566) || (fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) < 1.4442)) BE =true;
	    if(fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) > 1.566) EE =true;

	    tree->Fill();

	  } // pt
	} // trig matching
      } // good electrons
    } // event Loop

    file[jentry]->Write();
    file[jentry]->Close();
    printf ("sum of weights: %f \n",sum_weights);
  } // file Loop
}

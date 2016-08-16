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

void data13TeV() {
  //TFile f1("/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/nTuples/analysis_4March/SE_Run2015.root")
  TFile f1("eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/SE_2015.root");
  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  vector<float>   *genPostFSR_Pt;
  vector<float>   *genPostFSR_Eta;
  vector<float>   *genPostFSR_Rap;
  vector<float>   *genPostFSR_Phi;
  vector<float>   *genPostFSR_En;
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
  Int_t           nPV;
  vector<double>  *pt_Ele23;
  vector<double>  *eta_Ele23;
  vector<double>  *phi_Ele23;

  genPostFSR_Pt = 0;
  genPostFSR_Eta = 0;
  genPostFSR_Rap = 0;
  genPostFSR_Phi = 0;
  genPostFSR_En = 0;
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
  pt_Ele23 = 0;
  eta_Ele23 = 0;
  phi_Ele23 = 0;

  T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
  T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
  T1->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
  T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
  T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
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
  T1->SetBranchAddress("nPV", &nPV);
  T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
  T1->SetBranchAddress("pt_Ele23", &pt_Ele23);
  T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
  T1->SetBranchAddress("phi_Ele23", &phi_Ele23);

  TFile *file = new TFile("MediumId_Calibrated/SingleElectron_UnScaleCorr.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  int mediumId, dzMaskId, passId;
  bool trigMatch_lead, trigMatch_slead;
  double dR, dR1, dR2;
  TLorentzVector ele1,ele2,dielectron;
  TLorentzVector recoElec;

  // Branch variable declaration
  double primVtx;
  double Ele1PT, Ele2PT, Ele1Eta, Ele2Eta, Ele1Phi, Ele2Phi, Ele1Enr, Ele2Enr;
  double ZMass, ZPT, ZEta, ZRapidity, ZPhi;
  bool BB, BE, EE;

  // Branch declaration
  tree->Branch("primVtx", &primVtx, "primVtx/D");
  tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
  tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
  tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
  tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
  tree->Branch("Ele1Phi", &Ele1Phi, "Ele1Phi/D");
  tree->Branch("Ele2Phi", &Ele2Phi, "Ele2Phi/D");
   tree->Branch("Ele1Enr", &Ele1Enr, "Ele1Enr/D");
    tree->Branch("Ele2Enr", &Ele2Enr, "Ele2Enr/D");
  tree->Branch("ZMass", &ZMass, "ZMass/D");
  tree->Branch("ZPT", &ZPT, "ZPT/D");
  tree->Branch("ZEta", &ZEta, "ZEta/D");
  tree->Branch("ZRapidity", &ZRapidity, "ZRapidity/D");
  tree->Branch("ZPhi", &ZPhi, "ZPhi/D");
  tree->Branch("BB", &BB, "BB/B");
  tree->Branch("BE", &BE, "BE/B");
  tree->Branch("EE", &EE, "EE/B");

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi;vector <double> newscEta;
  vector <double> elePt; vector <double> eleEta; vector <double> eleEnr; vector <double> elePhi; vector <double> elescEta;
  TClonesArray *eleP4 = new TClonesArray("TLorentzVector", 1500);

  int nentries = T1->GetEntries();
  //int nentries = 10000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    T1->GetEntry(jentry);

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    trigMatch_lead = false; trigMatch_slead = false;
    dR = 0.; dR1 = 0.; dR2 = 0.;
    BB=false; BE=false; EE=false;
    mediumId = 0; dzMaskId = 0; passId =0;

    newelePt.clear(); neweleEta.clear(); newelePhi.clear(); neweleEnr.clear();newscEta.clear();
    elePt.clear(); eleEta.clear(); elePhi.clear(); eleEnr.clear();elescEta.clear();
    recoElec.SetPtEtaPhiE(0.,0.,0.,0.);
    eleP4->Clear();

    primVtx = nPV;

    if(!Ele23_WPLoose) continue;

    for(int i=0;i<ptElec->size();i++){

      recoElec.SetPtEtaPhiE(ptElec->at(i),etaElec->at(i), phiElec->at(i), energyElec->at(i));
      new ((*eleP4)[i]) TLorentzVector(recoElec);
      TLorentzVector* fourmom = (TLorentzVector*) eleP4->At(i);

   //  if(fourmom->Eta() < 1.4442) *fourmom *= 0.9994;
   //   if(fourmom->Eta() > 1.566) *fourmom *= 1.0034;

      if(fourmom->Eta() < 1.4442) *fourmom *= 1.0;
      if(fourmom->Eta() > 1.566) *fourmom *= 1.0;

     //   mediumId = passHEEPId->at(i);
       mediumId = passMediumId->at(i);

       passId = mediumId;

//      if (ptElec->at(i) < 500.) { passId = mediumId; }
//     else { passId = dzMaskId; }

      if(passId) {
	if(fabs(etaSC->at(i)) < 2.5 && !(fabs(etaSC->at(i)) > 1.4442 && fabs(etaSC->at(i)) < 1.566)){

	  newelePt.push_back(fourmom->Pt());
	  neweleEta.push_back(fourmom->Eta());
	  neweleEnr.push_back(fourmom->Energy());
	  newelePhi.push_back(fourmom->Phi());
          newscEta.push_back(etaSC->at(i));
	}
      }
    }

    // Sorting
    int index[newelePt.size()];
    float pt[newelePt.size()];

    for(unsigned int el=0; el<newelePt.size(); el++)
    { 
      pt[el]=newelePt.at(el);
    }

    int size = sizeof(pt)/sizeof(pt[0]);
    TMath::Sort(size,pt,index,true);

    if(newelePt.size()==2){

      for(unsigned int j = 0; j < pt_Ele23->size(); j++){

	double dR1_comp = 1000.;
	double dR2_comp = 1000.;

	dR1 = deltaR(neweleEta.at(index[0]), newelePhi.at(index[0]), eta_Ele23->at(j), phi_Ele23->at(j));
	dR2 = deltaR(neweleEta.at(index[1]), newelePhi.at(index[1]), eta_Ele23->at(j), phi_Ele23->at(j));

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

	if(newelePt.at(index[0]) > 30. && newelePt.at(index[1]) > 10.) {

	  ele1.SetPtEtaPhiE(newelePt.at(index[0]),neweleEta.at(index[0]),newelePhi.at(index[0]),neweleEnr.at(index[0]));
	  ele2.SetPtEtaPhiE(newelePt.at(index[1]),neweleEta.at(index[1]),newelePhi.at(index[1]),neweleEnr.at(index[1]));
	  dielectron=ele1+ele2;

	  Ele1PT  = newelePt.at(index[0]);
	  Ele2PT  = newelePt.at(index[1]);
	  Ele1Eta = neweleEta.at(index[0]);
	  Ele2Eta = neweleEta.at(index[1]);
	  Ele1Phi = newelePhi.at(index[0]);
	  Ele2Phi = newelePhi.at(index[1]);

	  ZMass = dielectron.M();
	  ZPT = dielectron.Pt();
	  ZEta = dielectron.Eta();
	  ZRapidity = dielectron.Rapidity();
	  ZPhi = dielectron.Phi();
/*
	  if(fabs(neweleEta.at(index[0])) < 1.4442 && fabs(neweleEta.at(index[1])) < 1.4442) BB = true;
	  if((fabs(neweleEta.at(index[0])) < 1.4442 && fabs(neweleEta.at(index[1])) > 1.566) || (fabs(neweleEta.at(index[0])) > 1.566 && fabs(neweleEta.at(index[1])) < 1.4442)) BE =true;
	  if(fabs(neweleEta.at(index[0])) > 1.566 && fabs(neweleEta.at(index[1])) > 1.566) EE =true;
*/

          if(fabs(newscEta.at(index[0])) < 1.4442 && fabs(newscEta.at(index[1])) < 1.4442) BB = true;
          if((fabs(newscEta.at(index[0])) < 1.4442 && fabs(newscEta.at(index[1])) > 1.566) || (fabs(newscEta.at(index[0])) > 1.566 && fabs(newscEta.at(index[1])) < 1.4442)) BE =true;
          if(fabs(newscEta.at(index[0])) > 1.566 && fabs(newscEta.at(index[1])) > 1.566) EE =true;
	  tree->Fill();

	} // pt
      } // trig matching
    } // size==2
  } // event

  file->Write();
  file->Close();
}

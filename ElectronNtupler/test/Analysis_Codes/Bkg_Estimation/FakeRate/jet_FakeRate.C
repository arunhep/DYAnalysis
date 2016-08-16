#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void jet_FakeRate() {

  TFile f1("/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/Photon_2015.root");
  //TFile f1("/tmp/rchawla/Photon_2015.root");
  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  vector<double>  *etSPhoHLT;
  vector<double>  *etaSPhoHLT;
  vector<double>  *phiSPhoHLT;
  Int_t           singlePhoton;
  Int_t           prescalePhoton;
  vector<int>     *expectedMissingInnerHits;
  vector<double>  *metPt;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *etaSC;
  vector<int>     *passMediumId;

  etSPhoHLT = 0;
  etaSPhoHLT = 0;
  phiSPhoHLT = 0;
  expectedMissingInnerHits = 0;
  metPt = 0;
  ptElec = 0;
  etaElec = 0;
  phiElec = 0;
  energyElec = 0;
  etaSC = 0;
  passMediumId = 0; 

  T1->SetBranchStatus("*", 0);
  T1->SetBranchStatus("etSPhoHLT", 1);
  T1->SetBranchStatus("etaSPhoHLT", 1);
  T1->SetBranchStatus("phiSPhoHLT", 1);
  T1->SetBranchStatus("singlePhoton", 1);
  T1->SetBranchStatus("prescalePhoton", 1);
  T1->SetBranchStatus("expectedMissingInnerHits", 1);
  T1->SetBranchStatus("metPt", 1);
  T1->SetBranchStatus("ptElec", 1);
  T1->SetBranchStatus("etaElec", 1);
  T1->SetBranchStatus("phiElec", 1);
  T1->SetBranchStatus("energyElec", 1);
  T1->SetBranchStatus("etaSC", 1);
  T1->SetBranchStatus("passMediumId", 1);

  T1->SetBranchAddress("etSPhoHLT", &etSPhoHLT); 
  T1->SetBranchAddress("etaSPhoHLT", &etaSPhoHLT); 
  T1->SetBranchAddress("phiSPhoHLT", &phiSPhoHLT); 
  T1->SetBranchAddress("singlePhoton", &singlePhoton);
  T1->SetBranchAddress("prescalePhoton", &prescalePhoton);
  T1->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits);
  T1->SetBranchAddress("metPt", &metPt);
  T1->SetBranchAddress("ptElec", &ptElec);
  T1->SetBranchAddress("etaElec", &etaElec);
  T1->SetBranchAddress("phiElec", &phiElec);
  T1->SetBranchAddress("energyElec", &energyElec);
  T1->SetBranchAddress("etaSC", &etaSC);
  T1->SetBranchAddress("passMediumId", &passMediumId);

  //TFile *file = new TFile("fakeRate_Photon_noMET_noMissHits.root", "recreate");
  TFile *file = new TFile("fakeRate_Photon_check.root", "recreate");

  int count;
  bool passKin;
  int mediumId, passId;
  double dR;
  vector <double> newelePt; vector <double> neweleEta;  vector <double> newscEta; vector <double> neweleMediumId;

  //Double_t x1bin[10] = {10,15,20,25,30,40,50,70,100,10000};
  Double_t x1bin[15] = {10,15,20,25,30,40,50,70,100,150,200,300,400,500,10000};
  int nbins = 14;

  TH1F *et_HLT = new TH1F("et_HLT", "et_HLT", 100, 0, 700);
  TH1F *et_HLT_preScale = new TH1F("et_HLT_preScale", "et_HLT_preScale", 100, 0, 700);

  TH1F *numPt      = new TH1F("numPt", "numPt", nbins, x1bin);
  TH1F *numPt_BRL  = new TH1F("numPt_BRL", "numPt_BRL", nbins, x1bin);
  TH1F *numPt_ECAP = new TH1F("numPt_ECAP", "numPt_ECAP", nbins, x1bin);

  TH1F *denPt      = new TH1F("denPt", "denPt", nbins, x1bin);
  TH1F *denPt_BRL  = new TH1F("denPt_BRL", "denPt_BRL", nbins, x1bin);
  TH1F *denPt_ECAP = new TH1F("denPt_ECAP", "denPt_ECAP", nbins, x1bin);

  et_HLT->Sumw2(); et_HLT_preScale->Sumw2();
  numPt->Sumw2(); denPt->Sumw2();
  numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
  numPt_ECAP->Sumw2(); denPt_ECAP->Sumw2();

  int nentries = T1->GetEntries();
  //int nentries = 50000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    T1->GetEntry(jentry);

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    // Sorting
    int index[ptElec->size()];
    float pt[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++) {
      pt[el]=ptElec->at(el); }

    int size = sizeof(pt)/sizeof(pt[0]);
    TMath::Sort(size,pt,index,true);

    count = 0;
    passKin = false;
    mediumId = 0; passId = 0;
    dR = 0.;
    newelePt.clear(); neweleEta.clear(); newscEta.clear(); neweleMediumId.clear();

    //if(metPt->at(0) > 10.) continue;
    if(!singlePhoton) continue;

    et_HLT->Fill(etSPhoHLT->at(0));
    et_HLT_preScale->Fill(etSPhoHLT->at(0),prescalePhoton);

    for(int i=0;i<ptElec->size();i++){

      mediumId = passMediumId->at(index[i]);
      passKin = (fabs(etaSC->at(index[i])) < 2.5 && !(fabs(etaSC->at(index[i])) > 1.4442 && fabs(etaSC->at(index[i])) < 1.566));

      if(mediumId && passKin) count++;
    }

    if(count <= 1){
      for(unsigned int j=0;j<ptElec->size();j++){
	//if(expectedMissingInnerHits->at(index[j]) == 0){

	newelePt.push_back(ptElec->at(index[j]));
	neweleEta.push_back(etaElec->at(index[j]));
	newscEta.push_back(etaSC->at(index[j]));
	neweleMediumId.push_back(passMediumId->at(index[j]));
	//}
      }
    }

    for(unsigned int k=0;k<newelePt.size();k++){
      if(newelePt.at(k)  > 10000.) continue;

      if(newelePt.at(k) > 10000.) cout<<"Check"<<endl;

      denPt->Fill(newelePt.at(k),prescalePhoton);

      if(fabs(newscEta.at(k)) < 1.4442) denPt_BRL->Fill(newelePt.at(k),prescalePhoton);
      else if(fabs(newscEta.at(k)) > 1.566 && fabs(newscEta.at(k)) < 2.5) denPt_ECAP->Fill(newelePt.at(k),prescalePhoton);

      passId = neweleMediumId.at(k);
      if(passId){

	numPt->Fill(newelePt.at(k),prescalePhoton);

	if(fabs(newscEta.at(k)) < 1.4442) numPt_BRL->Fill(newelePt.at(k),prescalePhoton);
	else if(fabs(newscEta.at(k)) > 1.566 && fabs(newscEta.at(k)) < 2.5) numPt_ECAP->Fill(newelePt.at(k),prescalePhoton);

      } // ID
    }// elePt size
  } // event

  file->Write();
  file->Close();
}

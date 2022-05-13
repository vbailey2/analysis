float subtractTower(float towerE, float towerUE, float towerPhi, float v2, float psi2, int doFlow = 1){
  
  if(doFlow){
    towerUE += towerUE * (2*v2*TMath::Cos(2*(towerPhi - psi2)));
  }
  
  return towerE-towerUE;

}

int getEtaBin(float eta){
  for(int i = 0; i < 24; i++){
    if(eta < -1.1+(i+1)*2.2/24) return i;
  }
  return -1;
}

int getPhiBin(float phi){
  for(int i = 0; i < 64; i++){
    if(phi < (i+1)*2*TMath::Pi()/64) return i;
  }
  return -1;
}

void getFluctuations(string in = "/sphenix/user/vbailey/MDC2HIJINGFluctuationNtuples/output.root", string out ="testhist", int subtraction = 0, int nEvents = 100000){
  //subtraction: 0= no subtraction, 1= subtraction without flow, 2= subtraction with flow
  int doFlow = (subtraction == 2) ? 1 : 0;
  
  //number of different window sizes we are looking at
  int nsizes = 3; //window sizes: 1x1, 4x3 (R = 0.2 jet), 7x7 (R = 0.4 jet)

  //read in file

  TChain* chain = new TChain("tree");
  chain->Add(in.c_str());
  if (chain == 0) return;
  Long64_t nentries = (nEvents > 0) ? nEvents+1 : chain->GetEntries();
  if(nentries == 0){
    std::cout<<"No entries"<<std::endl;
    return;
  }

  int m_event;
  std::vector<float> *m_towerCEMC_e = 0;
  std::vector<float> *m_towerCEMC_eta = 0;
  std::vector<float> *m_towerCEMC_phi = 0;
  std::vector<int> *m_towerCEMC_etaBin = 0;
  std::vector<int> *m_towerCEMC_phiBin = 0;

  std::vector<float> *m_towerIHCAL_e = 0;
  std::vector<float> *m_towerIHCAL_eta = 0;
  std::vector<float> *m_towerIHCAL_phi = 0;
  std::vector<int> *m_towerIHCAL_etaBin = 0;
  std::vector<int> *m_towerIHCAL_phiBin = 0;

  std::vector<float> *m_towerOHCAL_e = 0;
  std::vector<float> *m_towerOHCAL_eta = 0;
  std::vector<float> *m_towerOHCAL_phi = 0;
  std::vector<int> *m_towerOHCAL_etaBin = 0;
  std::vector<int> *m_towerOHCAL_phiBin = 0;

  std::vector<float> *m_UE0 = 0;
  std::vector<float> *m_UE1 = 0;
  std::vector<float> *m_UE2 = 0;

  float totalCEMC;
  float totalIHCAL;
  float totalOHCAL;

  float m_v2;
  double m_b;
  float m_psi2;

  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *pt = 0;
  vector<float> *e = 0;
  vector<float> *truthEta = 0;
  vector<float> *truthPhi = 0;
  vector<float> *truthPt = 0;
  vector<float> *truthE = 0;

  chain->SetBranchAddress("m_event", &m_event);
  chain->SetBranchAddress("m_towerCEMC_e", &m_towerCEMC_e);
  chain->SetBranchAddress("m_towerCEMC_eta", &m_towerCEMC_eta);
  chain->SetBranchAddress("m_towerCEMC_phi", &m_towerCEMC_phi);
  chain->SetBranchAddress("m_towerCEMC_etaBin", &m_towerCEMC_etaBin);
  chain->SetBranchAddress("m_towerCEMC_phiBin", &m_towerCEMC_phiBin);

  chain->SetBranchAddress("m_towerIHCAL_e", &m_towerIHCAL_e);
  chain->SetBranchAddress("m_towerIHCAL_eta", &m_towerIHCAL_eta);
  chain->SetBranchAddress("m_towerIHCAL_phi", &m_towerIHCAL_phi);
  chain->SetBranchAddress("m_towerIHCAL_etaBin", &m_towerIHCAL_etaBin);
  chain->SetBranchAddress("m_towerIHCAL_phiBin", &m_towerIHCAL_phiBin);

  chain->SetBranchAddress("m_towerOHCAL_e", &m_towerOHCAL_e);
  chain->SetBranchAddress("m_towerOHCAL_eta", &m_towerOHCAL_eta);
  chain->SetBranchAddress("m_towerOHCAL_phi", &m_towerOHCAL_phi);
  chain->SetBranchAddress("m_towerOHCAL_etaBin", &m_towerOHCAL_etaBin);
  chain->SetBranchAddress("m_towerOHCAL_phiBin", &m_towerOHCAL_phiBin);

  chain->SetBranchAddress("m_UE0", &m_UE0);
  chain->SetBranchAddress("m_UE1", &m_UE1);
  chain->SetBranchAddress("m_UE2", &m_UE2);

  chain->SetBranchAddress("m_totalCEMC",&totalCEMC);
  chain->SetBranchAddress("m_totalIHCAL",&totalIHCAL);
  chain->SetBranchAddress("m_totalOHCAL",&totalOHCAL);

  chain->SetBranchAddress("m_v2", &m_v2);
  chain->SetBranchAddress("m_b", &m_b);
  chain->SetBranchAddress("m_psi2", &m_psi2);

  chain->SetBranchAddress("eta",&eta);
  chain->SetBranchAddress("phi",&phi);
  chain->SetBranchAddress("pt",&pt);
  chain->SetBranchAddress("e",&e);
  chain->SetBranchAddress("truthEta",&truthEta);
  chain->SetBranchAddress("truthPhi",&truthPhi);
  chain->SetBranchAddress("truthPt",&truthPt);
  chain->SetBranchAddress("truthE",&truthE);


  //create histograms
  TH1D *h_Et = new TH1D("h_Et","",200,-10,10);
  TH2D *h_AvgE[nsizes][3];//average energy for a given window size: n window sizes, 3 subdetectors
  TH2D *h_STD[nsizes][3];//STD of energy for a given window size: n window sizes, 3 subdetectors  
  TH2D *h_AvgEtot[nsizes];
  TH2D *h_STDtot[nsizes];
  TH3D *h_etaPhiEt[nsizes][3];

  TH2D *h_AvgE_TE[nsizes][3];
  TH2D *h_AvgE_IP[nsizes][3];
  TH2D *h_STD_TE[nsizes][3];
  TH2D *h_STD_IP[nsizes][3];

  TH2D *h_AvgEtot_TE[nsizes];
  TH2D *h_STDtot_TE[nsizes];
  TH2D *h_AvgEtot_IP[nsizes];
  TH2D *h_STDtot_IP[nsizes];
  
  TH2D *h_window1x1Etot[3];
  TH2D *h_window3x4Etot[3];
  TH2D *h_window7x7Etot[3];


  double meanEt1x1sub[3];
  double meanEt3x4sub[3];
  double meanEt7x7sub[3];
  double STDEt1x1sub[3];
  double STDEt3x4sub[3];
  double STDEt7x7sub[3];
  double meanEt1x1;
  double meanEt3x4;
  double meanEt7x7;
  double STDEt1x1;
  double STDEt3x4;
  double STDEt7x7;
  TH2D *h_EventEt[3];

  string sizes[] = {"1x1","3x4","7x7"};
  string detectors[] = {"EMCal","IHCal","OHCal"};
  for(int i = 0; i < nsizes; i++){
    for(int j = 0; j < 3; j++){
      h_AvgE_IP[i][j]  = new TH2D(Form("h_AvgE%i_layer%i_IP",i,j),Form("Mean Window Energy vs Impact Parameter. Size %s Subdetector %s",sizes[i].c_str(),detectors[j].c_str()),200,0,20,8000,-10,40);
      h_STD_IP[i][j]  = new TH2D(Form("h_STD%i_layer%i_IP",i,j),Form("Standard Deviation Window Energy vs Impact Parameter. Size %s Subdetector %s",sizes[i].c_str(),detectors[j].c_str()),200,0,20,1400,0,7);
      h_AvgE_TE[i][j]  = new TH2D(Form("h_AvgE%i_layer%i_TE",i,j),Form("Mean Window Energy vs Total Energy. Size %s Subdetector %s",sizes[i].c_str(),detectors[j].c_str()),200,0,2500,8000,-10,40);
      h_STD_TE[i][j]  = new TH2D(Form("h_STD%i_layer%i_TE",i,j),Form("Standard Deviation Window Energy vs Total Energy. Size %s Subdetector %s",sizes[i].c_str(),detectors[j].c_str()),200,0,2500,1400,0,7);
    }
    h_AvgEtot_IP[i] = new TH2D(Form("h_AvgEtot%i_IP",i),Form("Mean Winodw Energy vs Impact Parameter. Window Size %s",sizes[i].c_str()),200,0,20,200,-1,1);
    h_STDtot_IP[i] = new TH2D(Form("h_STDtot%i_IP",i),Form("Standard Deviation Winodw Energy vs Impact Parameter. Window Size %s",sizes[i].c_str()),200,0,20,1400,0,7);
    h_AvgEtot_TE[i] = new TH2D(Form("h_AvgEtot%i_TE",i),Form("Mean Winodw Energy vs Total Energy. Window Size %s",sizes[i].c_str()),200,0,2500,200,-1,1);
    h_STDtot_TE[i] = new TH2D(Form("h_STDtot%i_TE",i),Form("Standard Deviation Winodw Energy vs Total Energy. Window Size %s",sizes[i].c_str()),200,0,2500,1400,0,7);
  }
  for(int isys = 0; isys < 3; isys++){
    h_window1x1Etot[isys] = new TH2D(Form("h_window1x1Etot_sub%i",isys),Form("1x1 Window Energy Subdetector %i",isys),24,-1.1,1.1,64,0,2*TMath::Pi());
    h_window3x4Etot[isys] = new TH2D(Form("h_window3x4Etot_sub%i",isys),Form("3x4 Window Energy Subdetector %i",isys),8,-1.1,1.1,16,0,2*TMath::Pi());
    h_window7x7Etot[isys] = new TH2D(Form("h_window7x7Etot_sub%i",isys),Form("7x7 Window Energy Subdetector %i",isys),3,-1.1+2.2/24,1.1-4.4/24,9,0,(2-1/64)*TMath::Pi());
  }
  //hists for total calo
  h_EventEt[0] = new TH2D("h_EventEt_1x1","",24,-1.1,1.1,64,0,2*TMath::Pi());
  h_EventEt[1] = new TH2D("h_EventEt_3x4","",8,-1.1,1.1,16,0,2*TMath::Pi());
  h_EventEt[2] = new TH2D("h_EventEt_7x7","",3,-1.1+2.2/24,1.1-4.4/24,9,0,(2-1/64)*TMath::Pi());

  //loop through all the entries in our file and fill the histograms

  std::cout<<"running "<<nentries<<" events"<<std::endl;
  for(Long64_t jentry=1; jentry<nentries;jentry++){                 int e=jentry;
    if(jentry%1000 == 0) std::cout<<"Entry: "<<jentry<<std::endl;
    
    chain->GetEntry(jentry);

    //get the total energy in all three subdetectors
    float totalEnergy = totalCEMC+totalIHCAL+totalOHCAL;
    float towerEt[64][24] = {}; //there are 64 towers in phi x 24 towers in eta
    
    double meantot1x1[3]={0,0,0}, meantot3x4[3]={0,0,0}, meantot7x7[3]={0,0,0};
    double rmstot1x1[3]={0,0,0}, rmstot3x4[3]={0,0,0}, rmstot7x7[3]={0,0,0};
    //loop through each subdetector: 0: EMCal, 1: inner HCal, 2: outer HCal
    for(int isys = 0; isys < 3; isys++){
        
      //get the number of towers in an event for the given subdetector 
      int nTower = 0;
      if(isys == 0) nTower = m_towerCEMC_e->size();
      else if(isys == 1) nTower = m_towerIHCAL_e->size();
      else if(isys == 2) nTower = m_towerOHCAL_e->size();

      //get the energy for each tower (and subtract the underlying event if using subtraction)
      for(int i = 0; i < nTower; i++){
	float eta;
	float phi;
	int etabin;
	int phibin;
	float energy;
	if(isys == 0){
	  eta = m_towerCEMC_eta->at(i);
	  phi = m_towerCEMC_phi->at(i);
	  energy = (subtraction == 0) ? m_towerCEMC_e->at(i) : subtractTower(m_towerCEMC_e->at(i), m_UE0->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
	}
	else if(isys== 1){
          eta = m_towerIHCAL_eta->at(i);
          phi = m_towerIHCAL_phi->at(i);
	  energy = (subtraction == 0) ? m_towerIHCAL_e->at(i) : subtractTower(m_towerIHCAL_e->at(i), m_UE1->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
	}
	if(isys== 2){
          eta = m_towerOHCAL_eta->at(i);
          phi = m_towerOHCAL_phi->at(i);
	  energy = (subtraction == 0) ? m_towerOHCAL_e->at(i) : subtractTower(m_towerOHCAL_e->at(i), m_UE2->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
	}
	etabin = getEtaBin(eta);
	phibin = getPhiBin(phi);
	towerEt[phibin][etabin] = energy/cosh(eta);
	h_Et->Fill(towerEt[phibin][etabin]);
	h_window1x1Etot[isys]->Fill(eta, phi, towerEt[phibin][etabin]);
	h_window3x4Etot[isys]->Fill(eta, phi, towerEt[phibin][etabin]);
	h_window7x7Etot[isys]->Fill(eta, phi, towerEt[phibin][etabin]);
      }

    //calculate mean and RMS of window Et
    double rrmss;
    for (int i=1;i<25;i++){
        for (int j=1;j<65;j++){
            meantot1x1[isys] += h_window1x1Etot[isys]->GetBinContent(i,j);
            rrmss = h_window1x1Etot[isys]->GetBinContent(i,j);
            rmstot1x1[isys] += rrmss*rrmss;
        }
    }
    for (int i=1;i<9;i++){
        for (int j=1;j<17;j++){
            meantot3x4[isys] += h_window3x4Etot[isys]->GetBinContent(i,j);
            rrmss = h_window3x4Etot[isys]->GetBinContent(i,j);
            rmstot3x4[isys] += rrmss*rrmss;
        }
    }
    
    for (int i=1;i<4;i++){
      for (int j=1;j<10;j++){
	meantot7x7[isys] += h_window7x7Etot[isys]->GetBinContent(i,j); 
	rrmss = h_window7x7Etot[isys]->GetBinContent(i,j); 
	rmstot7x7[isys] += rrmss*rrmss;
      }
    }
    meanEt1x1sub[isys] = meantot1x1[isys] / 1536;
    meanEt3x4sub[isys] = meantot3x4[isys] / 128;
    meanEt7x7sub[isys] = meantot7x7[isys] / 27;
    STDEt1x1sub[isys] = (rmstot1x1[isys] / 1536);
    STDEt3x4sub[isys] = (rmstot3x4[isys] / 128);
    STDEt7x7sub[isys] = (rmstot7x7[isys] / 27);

    h_AvgE_IP[0][isys]->Fill(m_b,meanEt1x1sub[isys]);
    h_AvgE_IP[1][isys]->Fill(m_b,meanEt3x4sub[isys]);
    h_AvgE_IP[2][isys]->Fill(m_b,meanEt7x7sub[isys]);
    h_STD_IP[0][isys]->Fill(m_b,TMath::Sqrt(STDEt1x1sub[isys] - meanEt1x1sub[isys]*meanEt1x1sub[isys] ) );
    h_STD_IP[1][isys]->Fill(m_b,TMath::Sqrt(STDEt3x4sub[isys] - meanEt3x4sub[isys]*meanEt3x4sub[isys] ) );
    h_STD_IP[2][isys]->Fill(m_b,TMath::Sqrt(STDEt7x7sub[isys] - meanEt7x7sub[isys]*meanEt7x7sub[isys] ) );

    h_AvgE_TE[0][isys]->Fill(totalEnergy,meanEt1x1sub[isys]);
    h_AvgE_TE[1][isys]->Fill(totalEnergy,meanEt3x4sub[isys]);
    h_AvgE_TE[2][isys]->Fill(totalEnergy,meanEt7x7sub[isys]);
    h_STD_TE[0][isys]->Fill(totalEnergy,TMath::Sqrt(STDEt1x1sub[isys] - meanEt1x1sub[isys]*meanEt1x1sub[isys] ) );
    h_STD_TE[1][isys]->Fill(totalEnergy,TMath::Sqrt(STDEt3x4sub[isys] - meanEt3x4sub[isys]*meanEt3x4sub[isys] ) );
    h_STD_TE[2][isys]->Fill(totalEnergy,TMath::Sqrt(STDEt7x7sub[isys] - meanEt7x7sub[isys]*meanEt7x7sub[isys] ) );
    }


    for (int i=0; i<3;i++){
        meanEt1x1 += meantot1x1[i];
        meanEt3x4 += meantot3x4[i];
        meanEt7x7 += meantot7x7[i];
    }
    meanEt1x1 /= 1536;
    meanEt3x4 /= 128;
    meanEt7x7 /= 27;
    double rrmss=0;
    double rms2tot1x1=0; double rms2tot3x4=0; double rms2tot7x7=0;
    for (int i=1;i<25;i++){
        for (int j=1;j<65;j++){
            for (int isys=0;isys<3;isys++){
                rrmss += h_window1x1Etot[isys]->GetBinContent(i,j);
            }
            h_EventEt[0]->Fill(i,j,rrmss);
            rms2tot1x1 += rrmss*rrmss;
            rrmss=0;
        }
    }
    rrmss=0;
    for (int i=1;i<9;i++){
        for (int j=1;j<17;j++){
            for (int isys=0;isys<3;isys++){
                rrmss += h_window3x4Etot[isys]->GetBinContent(i,j);
            }
            h_EventEt[1]->Fill(i,j,rrmss);
            rms2tot3x4 += rrmss*rrmss;
            rrmss=0;
        }
    }
    rrmss=0;
    for (int i=1;i<4;i++){
        for (int j=1;j<10;j++){
            for (int isys=0;isys<3;isys++){
                rrmss += h_window7x7Etot[isys]->GetBinContent(i,j); 
            }
            h_EventEt[2]->Fill(i,j,rrmss);
            rms2tot7x7 += rrmss*rrmss;
            rrmss=0;
        }
    }
    STDEt1x1 = (rms2tot1x1 / 1536);
    STDEt3x4 = (rms2tot3x4 / 128);
    STDEt7x7 = (rms2tot7x7 / 27);
    
    h_AvgEtot_IP[0]->Fill(m_b,meanEt1x1 / 1536);
    h_AvgEtot_IP[1]->Fill(m_b,meanEt3x4 / 128);
    h_AvgEtot_IP[2]->Fill(m_b,meanEt7x7 / 27);
    h_STDtot_IP[0]->Fill(m_b, TMath::Sqrt(STDEt1x1 - meanEt1x1 * meanEt1x1) );
    h_STDtot_IP[1]->Fill(m_b, TMath::Sqrt(STDEt3x4 - meanEt3x4 * meanEt3x4) );
    h_STDtot_IP[2]->Fill(m_b, TMath::Sqrt(STDEt7x7 - meanEt7x7 * meanEt7x7) );

    h_AvgEtot_TE[0]->Fill(totalEnergy,meanEt1x1 / 1536);
    h_AvgEtot_TE[1]->Fill(totalEnergy,meanEt3x4 / 128);
    h_AvgEtot_TE[2]->Fill(totalEnergy,meanEt7x7 / 27);
    h_STDtot_TE[0]->Fill(totalEnergy, TMath::Sqrt(STDEt1x1 - meanEt1x1 * meanEt1x1) );
    h_STDtot_TE[1]->Fill(totalEnergy, TMath::Sqrt(STDEt3x4 - meanEt3x4 * meanEt3x4) );
    h_STDtot_TE[2]->Fill(totalEnergy, TMath::Sqrt(STDEt7x7 - meanEt7x7 * meanEt7x7) );



    for(int k = 0; k < 3; k++){
      h_window1x1Etot[k]->Reset();
      h_window3x4Etot[k]->Reset();
      h_window7x7Etot[k]->Reset();
      h_EventEt[k]->Reset();
    }
  }

  string s = "";
  if(subtraction == 0) s = "_unsubtracted";
  else if(subtraction == 1) s = "_subtracted_noflow";
  else if(subtraction == 2) s = "_subtracted_withflow";
  TFile *f = new TFile(Form("%s%s.root",out.c_str(),s.c_str()),"RECREATE");
 
  for(int i = 0; i < nsizes; i++){
    h_AvgEtot_IP[i]->Write();
    h_STDtot_IP[i]->Write();
    h_AvgEtot_TE[i]->Write();
    h_STDtot_TE[i]->Write();
    for(int j = 0; j < 3; j++){
      h_STD_IP[i][j]->Write();
      h_AvgE_IP[i][j]->Write();
      h_STD_TE[i][j]->Write();
      h_AvgE_TE[i][j]->Write();
    }
  }

  f->Close();
}

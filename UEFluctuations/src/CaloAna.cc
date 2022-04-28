#include "CaloAna.h"

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

//Jet background includes
#include <jetbackground/TowerBackground.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/EventHeaderv1.h>

#include <g4jets/JetMapv1.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <TLorentzVector.h>

using namespace std;

CaloAna::CaloAna(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
  , hm(nullptr)
  , outfile(nullptr)
  , m_T(nullptr)
{
}

CaloAna::~CaloAna()
{
  delete hm;
}

int CaloAna::Init(PHCompositeNode*)
{
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  m_T = new TTree("tree", "MyAnalysis Tree");
  m_T->Branch("m_event", &m_event);
  m_T->Branch("m_towerCEMC_e", &m_towerCEMC_e);
  m_T->Branch("m_towerCEMC_eta", &m_towerCEMC_eta);
  m_T->Branch("m_towerCEMC_phi", &m_towerCEMC_phi);
  m_T->Branch("m_towerCEMC_etaBin", &m_towerCEMC_etaBin);
  m_T->Branch("m_towerCEMC_phiBin", &m_towerCEMC_phiBin);

  m_T->Branch("m_towerIHCAL_e", &m_towerIHCAL_e);
  m_T->Branch("m_towerIHCAL_eta", &m_towerIHCAL_eta);
  m_T->Branch("m_towerIHCAL_phi", &m_towerIHCAL_phi);
  m_T->Branch("m_towerIHCAL_etaBin", &m_towerIHCAL_etaBin);
  m_T->Branch("m_towerIHCAL_phiBin", &m_towerIHCAL_phiBin);

  m_T->Branch("m_towerOHCAL_e", &m_towerOHCAL_e);
  m_T->Branch("m_towerOHCAL_eta", &m_towerOHCAL_eta);
  m_T->Branch("m_towerOHCAL_phi", &m_towerOHCAL_phi);
  m_T->Branch("m_towerOHCAL_etaBin", &m_towerOHCAL_etaBin);
  m_T->Branch("m_towerOHCAL_phiBin", &m_towerOHCAL_phiBin);

  m_T->Branch("m_totalCEMC", &m_totalCEMC);
  m_T->Branch("m_totalIHCAL", &m_totalIHCAL);
  m_T->Branch("m_totalOHCAL", &m_totalOHCAL);

  m_T->Branch("m_UE0", &m_UE0);
  m_T->Branch("m_UE1", &m_UE1);
  m_T->Branch("m_UE2", &m_UE2);

  m_T->Branch("m_v2", &m_v2);
  m_T->Branch("m_psi2", &m_psi2);
  //m_T->Branch("m_v2_true", &m_v2_true);
  //m_T->Branch("m_psi2_true", &m_psi2_true);

  m_T->Branch("m_b", &m_b);

  return 0;
}

int CaloAna::process_event(PHCompositeNode* topNode)
{
  // For the calorimeters we have the following node name convention
  // where detector is the calorimeter name (CEMC, HCALIN, HCALOUT)
  // this is the order in which they are reconstructed
  //  G4HIT_<detector>: G4 Hits
  //  G4CELL_<detector>: Cells (combined hits inside a cell - scintillator, eta/phi bin)
  //  TOWER_SIM_<detector>: simulated tower (before pedestal and noise)
  //  TOWER_RAW_<detector>: Raw Tower (adc/tdc values - from sims or real data)
  //  TOWER_CALIB_<detector>: Calibrated towers
  //  CLUSTER_<detector>: clusters

  process_towers(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_g4hits(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_g4cells(PHCompositeNode* topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  //std::cout<<"Starting event"<<std::endl;
  ostringstream nodename;
  ostringstream geonodename;

  std::string thisDetector = "CEMC";
 
  nodename.str("");
  nodename << "TOWER_CALIB_" << thisDetector;
  geonodename.str("");
  geonodename << "TOWERGEOM_" << thisDetector;
  //RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, geonodename.str().c_str());
  
  ///since EMCal is retowered we want to use the HCal geometry to get the 0.1x0.1 EMCal towers (to get the 0.025x0.025 towers use the commented out above)
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  
  if (!towergeom)
  {
    std::cout<<"No CEMC"<<std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //for 0.025x0.025 EMCal towers uncomment below
  //RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, nodename.str().c_str());
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
  if (towers)
  {
    m_totalCEMC = towers->getTotalEdep();
    // again pair of iterators to begin and end of tower map
    RawTowerContainer::ConstRange tower_range = towers->getTowers();
    for (RawTowerContainer::ConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; tower_iter++)

    {
      
      int phibin = tower_iter->second->get_binphi();
      int etabin = tower_iter->second->get_bineta();
      float phi = towergeom->get_phicenter(phibin);
      float eta = towergeom->get_etacenter(etabin);
      m_towerCEMC_e.push_back(tower_iter->second->get_energy());
      m_towerCEMC_eta.push_back(eta);
      m_towerCEMC_phi.push_back(phi);
      m_towerCEMC_etaBin.push_back(etabin);
      m_towerCEMC_phiBin.push_back(phibin);
    }
  }

  thisDetector = "HCALIN";
  nodename.str("");
  nodename << "TOWER_CALIB_" << thisDetector;
  geonodename.str("");
  geonodename << "TOWERGEOM_" << thisDetector;
  RawTowerGeomContainer* towergeomIH = findNode::getClass<RawTowerGeomContainer>(topNode, geonodename.str().c_str());
  if (!towergeomIH)
    {
      std::cout<<"No IHCAl"<<std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  RawTowerContainer* towersIH = findNode::getClass<RawTowerContainer>(topNode, nodename.str().c_str());
  if (towersIH)
    {
      m_totalIHCAL= towersIH->getTotalEdep();
      RawTowerContainer::ConstRange tower_range = towersIH->getTowers();
      for (RawTowerContainer::ConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; tower_iter++)

	{
	  int phibin = tower_iter->second->get_binphi();
	  int etabin = tower_iter->second->get_bineta();
	  float phi = towergeomIH->get_phicenter(phibin);
	  float eta = towergeomIH->get_etacenter(etabin);
	  m_towerIHCAL_e.push_back(tower_iter->second->get_energy());
	  m_towerIHCAL_eta.push_back(eta);
	  m_towerIHCAL_phi.push_back(phi);
	  m_towerIHCAL_etaBin.push_back(etabin);
          m_towerIHCAL_phiBin.push_back(phibin);
	}
    }

  thisDetector = "HCALOUT";
  nodename.str("");
  nodename << "TOWER_CALIB_" << thisDetector;
  geonodename.str("");
  geonodename << "TOWERGEOM_" << thisDetector;
  RawTowerGeomContainer* towergeomOH = findNode::getClass<RawTowerGeomContainer>(topNode, geonodename.str().c_str());
  if (!towergeomOH)
    {
      std::cout<<"No OHCAL"<<std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  RawTowerContainer* towersOH = findNode::getClass<RawTowerContainer>(topNode, nodename.str().c_str());
  if (towersOH)
    {
      m_totalOHCAL = towersOH->getTotalEdep();
      RawTowerContainer::ConstRange tower_range = towersOH->getTowers();
      for (RawTowerContainer::ConstIterator tower_iter = tower_range.first; tower_iter != tower_range.second; tower_iter++)

        {
          int phibin = tower_iter->second->get_binphi();
          int etabin = tower_iter->second->get_bineta();
          float phi = towergeomOH->get_phicenter(phibin);
          float eta = towergeomOH->get_etacenter(etabin);
          m_towerOHCAL_e.push_back(tower_iter->second->get_energy());
          m_towerOHCAL_eta.push_back(eta);
          m_towerOHCAL_phi.push_back(phi);
	  m_towerOHCAL_etaBin.push_back(etabin);
          m_towerOHCAL_phiBin.push_back(phibin);
	}
    }

  TowerBackground *background=nullptr;
  background = findNode::getClass<TowerBackground>(topNode, "TowerBackground_Sub2");
  if(!background){
    std::cout<<"Can't get background. Exiting"<<std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_UE0 = background->get_UE(0);
  m_UE1 = background->get_UE(1);
  m_UE2 = background->get_UE(2);

  m_v2 = background->get_v2();
  m_psi2 = background->get_Psi2();

  //get impact parameter
  EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
  m_b = event_header->get_floatval("bimp");


  m_T->Fill();
  
  m_towerCEMC_e.clear();
  m_towerCEMC_eta.clear();
  m_towerCEMC_phi.clear();
  m_towerCEMC_etaBin.clear();
  m_towerCEMC_phiBin.clear();
  m_towerIHCAL_e.clear();
  m_towerIHCAL_eta.clear();
  m_towerIHCAL_phi.clear();
  m_towerIHCAL_etaBin.clear();
  m_towerIHCAL_phiBin.clear();
  m_towerOHCAL_e.clear();
  m_towerOHCAL_eta.clear();
  m_towerOHCAL_phi.clear();
  m_towerOHCAL_etaBin.clear();
  m_towerOHCAL_phiBin.clear();
  m_UE0.clear();
  m_UE1.clear();
  m_UE2.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_clusters(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::End(PHCompositeNode* topNode)
{
  outfile->cd();
  m_T->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}

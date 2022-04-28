#ifndef CALOANA_H__
#define CALOANA_H__

#include <fun4all/SubsysReco.h>
#include <vector>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;

class CaloAna : public SubsysReco
{
 public:
  //! constructor
  CaloAna(const std::string &name = "CaloAna", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloAna();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }

 protected:
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;
  TFile *outfile;
  TTree *m_T;

  int m_event;
  std::vector<float> m_towerCEMC_e;
  std::vector<float> m_towerCEMC_eta;
  std::vector<float> m_towerCEMC_phi;
  std::vector<int> m_towerCEMC_etaBin;
  std::vector<int> m_towerCEMC_phiBin;

  std::vector<float> m_towerIHCAL_e;
  std::vector<float> m_towerIHCAL_eta;
  std::vector<float> m_towerIHCAL_phi;
  std::vector<int> m_towerIHCAL_etaBin;
  std::vector<int> m_towerIHCAL_phiBin;

  std::vector<float> m_towerOHCAL_e;
  std::vector<float> m_towerOHCAL_eta;
  std::vector<float> m_towerOHCAL_phi;
  std::vector<int> m_towerOHCAL_etaBin;
  std::vector<int> m_towerOHCAL_phiBin;


  float m_totalCEMC;
  float m_totalIHCAL;
  float m_totalOHCAL;

  std::vector<float> m_UE0;
  std::vector<float> m_UE1;
  std::vector<float> m_UE2;

  float m_v2;
  float m_psi2;

  float m_v2_true;
  float m_psi2_true;

  double m_b;

};

#endif

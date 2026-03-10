#ifndef __BUILDMETATOWERS_H__
#define __BUILDMETATOWERS_H__

//////////////////////////////////////////////////////////////////////////
//									//
//			Build Meta Towers				//
//									//
// 	Does what is says on the tin, take in calorimeters, applies a   //
// 	vertex based correction and combines the towers to make a	//
// 	single "shifted em+hadronic calorimeter" w/ HCAL bins		//
// 	Works with vandy skimmer data and Tower info			//
//									//
//	Authors: 	Ben Kimelman, Skaydi 				//
//	First commit: 	9 March 2026					//
//	This commit: 	10 March 2026					//
//	version: 	v1.0						//
//									//
//	Notes on this commit: First commit				//
//									//
//////////////////////////////////////////////////////////////////////////


//c++ classes

#include <math.h>
#include <vector>
#include <array>
#include <string>

//Vandy classes
#include <vandyclasses/Tower.h>
#include <vadyclasses/EventInfo.h>

//Calo classes 
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

//G4 objects
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>

//phool 
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/PHNodeOperation.h>

class BuildMetaTowers
{
	public:
		BuildMetaTowers(
				float radius = 1.245 /*average inner hcal radius*/, 
				std::string input="VandyClass" /*Vandy Class, TowerInfo or Array*/
			       ) 
		{
			//assume that if radius given is small it is given in m, else
			if(radius > 50 ) this->R = radius;
			if(radius < 50 ) this->R = radius * 100;
			this->T	= input;
		}
		~BuildMetaTowers(){};
		
		//struct for Array input 
		class TowerArrayEntry
		{
			double Energy;
			double phicenter;
			double etacenter;
			int index;
		};

		//Common Methods
		void RunMetaTowerBuilder(float z);
		//Fun4All load in 
		void GetEMCALTowers(
				TowerInfoContainerv2 EMCAL, 
				RawTowerGeomContainer_Cylinderv1 EMCAL_geom, 
				bool retower=false
				); //if using Fun4AllTower objects

		void GetHCALTowers(
				TowerInfoContainerv2 HCAL, 
				RawTowerGeomContainer_Cylinderv1 HCAL_geom, 
				bool outer=false
				); //if using Fun4AllTower objects
		//Vandy Class load in 
		void GetEMCALTowers(
				Tower EMCAL,
				bool retower=false
				);
		void GetHCALTowers(
				Tower HCAL,
				bool outer=false
				);
		//Just basic array 
		void GetEMCALTowers(
				std::vector<TowerArrayEntry> EMCAL,
			       	bool retower = false;
				);
		void GetHCALTowers(
				std::vector<TowerArrayEntry> HCAL,
				bool outer=false
				);
		//just returns configs
		void setRadius	( float r = 1.245 ) { this->R = r ; } ;
	       	float getRadius	( ) { return this->R; };
		
		void 		setInput
				( std::string input = "VandyClass" ) { this->T = input; }; 
		std::string 	getInput() { return this->T; }; 
		
		//return the final object
		std::array<TowerArrayEntry, 1536>* getMetaTowers() { return this->MetaTowers; };

	private:
		//Private variables
		float 		R { 1245 }; //IHCAL half radius in cm
		std::string 	T { "VandyClass" }; //What input Type to expect
		float 	R_OHCAL_Outer {269};
		float 	R_OHCAL_Inner {182};
		float 	R_IHCAL_Outer {137};
		float 	R_IHCAL_Inner {116};
		float 	R_EMCAL_Outer {116};
		float	R_EMCAL_Inner { 90};	
		std::array <TowerArrayEntry, 1536>* MetaTowers = new std::array <TowerArrayEntry, 1536> {}; 
		std::array <TowerArrayEntry, 1536>* EMReTowers = new std::array <TowerArrayEntry, 1536> {}; 
		std::array <TowerArrayEntry, 1536>* IHCaTowers = new std::array <TowerArrayEntry, 1536> {}; 
		std::array <TowerArrayEntry, 1536>* OHCaTowers = new std::array <TowerArrayEntry, 1536> {};
		
		//HCAL eta-phi physical geometry
		const std::array <double, 24> etaEdges 	
			{ -1.1000000, -1.0083333, -0.91666667, -0.82500000, -0.73333333, -0.64166667, 
				-0.55000000, -0.45833333, -0.36666667, -0.27500000, -0.18333333, -0.091666667, 
				1.1102230e-16, 0.091666667, 0.18333333, 0.27500000, 0.36666667, 0.45833333, 
				0.55000000, 0.64166667, 0.73333333, 0.82500000, 0.91666667, 1.0083333, 1.1000000 };

		const std::array <double, 64> phiEdges
			{ -0.053619781, 0.044554989, 0.14272976, 0.24090453, 0.33907930, 0.43725407, 
				0.53542884, 0.63360361, 0.73177838, 0.82995315, 0.92812792, 1.0263027, 
				1.1244775, 1.2226522, 1.3208270, 1.4190018, 1.5171765, 1.6153513, 1.7135261, 
				1.8117009, 1.9098756, 2.0080504, 2.1062252, 2.2043999, 2.3025747, 2.4007495, 
				2.4989242, 2.5970990, 2.6952738, 2.7934486, 2.8916233, 2.9897981, 3.0879729, 
				3.1861476, 3.2843224, 3.3824972, 3.4806720, 3.5788467, 3.6770215, 3.7751963, 
				3.8733710, 3.9715458, 4.0697206, 4.1678953, 4.2660701, 4.3642449, 4.4624197, 
				4.5605944, 4.6587692, 4.7569440, 4.8551187, 4.9532935, 5.0514683, 5.1496431, 
				5.2478178, 5.3459926, 5.4441674, 5.5423421, 5.6405169, 5.7386917, 5.8368664, 
				5.9350412, 6.0332160, 6.1313908, 6.2295655 };	
		
		std::array<double ,24> shiftedetaEdges {};

		//Private Methods
		double calculateEtaShift(double eta, double zVtx)
		{
			double z 	= R*std::sinh(eta);
			double zshift	= z-zVtx;
			zshift		= zshift / R;
			double etashift	= std::asinh(zshift);
			return etashift;
		}

		int calculateIndex(TowerArrayEntry tower) 
		{
			int index = tower.eta * 64 + tower.phi ;
			return index;
		}
		
		void decodeIndex(int index, int* ebr, int* pbr) 
		{
			
			int etabin	= index / 64;
			int phibin	= index % 64;
			ebr = &etabin;
			pbr = &phibin;
			return;
		}
		
		void shiftBinMins(double zVtx)
		{
			for(int i= 0; i<(int) etaEdges.size(); i++)
			{
				shiftedetaEdges[i]=calculateEtaShift(etaEdges[i], zVtx);
			}
			return;
		}
		void shiftEMCAL(TowerInfoContainer* emcal, bool retower);
		void shiftIHCAL(TowerInfoContainer* ihcal);
		void shiftOHCAL(TowerInfoContainer* ohcal);
		void shiftEMCAL(Tower emcal, bool retower);
		void shiftIHCAL(Tower ihcal);
		void shiftOHCAL(Tower ohcal);
		void shiftEMCAL(TowerArrayEntry emcal, bool retower);
		void shiftIHCAL(Tower ihcal);
		void shiftOHCAL(Tower ohcal);
		void addMetaTower(TowerEntryArray tower)
		{
			double eta_val 	= tower.eta;
			double phi_val	= tower.phi;
			int etaBin 	= findEtaBin(eta_val);
			int phiBin	= findPhiBin(phi_val);
			int N 		= calculateIndex(etaBin, phiBin);

			MetaTowers->at(N).E += tower.E;
			return;
		}
};
			
#endif

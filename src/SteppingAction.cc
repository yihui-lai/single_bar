#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DR_PMTSD.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
//#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();

using namespace std;
using namespace CLHEP;
/*
SteppingAction::SteppingAction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;

  config.readInto(core_material, "core_material"); 

  if (core_material == 0)
  {
	  config.readInto(toy_ly,	"toy_ly");
	  config.readInto(toy_decay,	"toy_decay");
	  config.readInto(toy_rise,	"toy_rise");
  }
  

}*/

int to_int(string name)
{
  int Result; // int which will contain the result
  stringstream convert(name);
  string dummy;
  convert >> dummy;
  convert >> Result;
  return Result;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::SteppingAction(DetectorConstruction *detectorConstruction,
                               const G4int &scint, const G4int &cher) : fDetectorConstruction(detectorConstruction),
                                                                        propagateScintillation(scint),
                                                                        propagateCerenkov(cher)
{
  maxtracklength = 500000. * mm;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::~SteppingAction()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void SteppingAction::UserSteppingAction(const G4Step *theStep)
{

  G4Track *theTrack = theStep->GetTrack();

  G4ParticleDefinition *particleType = theTrack->GetDefinition();

  G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
  G4StepPoint *thePostPoint = theStep->GetPostStepPoint();
  const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
  G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName = "";
  if (thePrePV)
    thePrePVName = thePrePV->GetName();
  G4String thePostPVName = "";
  if (thePostPV)
    thePostPVName = thePostPV->GetName();

  G4int nStep = theTrack->GetCurrentStepNumber();
  G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

  //        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------

  // get position
  G4double global_x = thePrePosition.x() / mm;
  G4double global_y = thePrePosition.y() / mm;
  G4double global_z = thePrePosition.z() / mm;

  G4double energy = theStep->GetTotalEnergyDeposit();
  G4double energyIon = energy - theStep->GetNonIonizingEnergyDeposit();
  G4double energyElec = 0;
  if (abs(TrPDGid) == 11) energyElec = energyIon;

  CreateTree::Instance()->depositedEnergyTotal += energy / GeV;
  CreateTree::Instance()->depositedIonEnergyTotal += energyIon / GeV;
  CreateTree::Instance()->depositedElecEnergyTotal += energyElec / GeV;

  bool outworld = ((theStep->GetPostStepPoint())->GetStepStatus()) == fWorldBoundary;
  if (outworld)
  {
    CreateTree::Instance()->depositedEnergyEscapeWorld += (theStep->GetPostStepPoint())->GetKineticEnergy() / GeV;
  }

  //------------- optical photon -------------
  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    //if optics
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
    float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / eV);
    //only consider 300 to 1000nm
    if (photWL > 1000 || photWL < 300)
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    else
    {
      //only consider Cerenkov and Scintillation
      if ((processName == "Cerenkov") || processName == "Scintillation")
      {
        //if reach the detector's window 
        if ((thePostPVName.contains("ecalDetWindowP_ff_")))
        {
          if (processName == "Scintillation")
          {
            //CreateTree::Instance()->h_photon_2D_receive_Scin->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_f_collect_Scin->Fill(photWL);
            CreateTree::Instance()->SDdetected_ff_S++;
          }
          if (processName == "Cerenkov")
          {
            //CreateTree::Instance()->h_photon_2D_receive_Ceren->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_f_collect_Ceren->Fill(photWL);
            CreateTree::Instance()->SDdetected_ff_C++;
          }
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        else if (thePostPVName.contains("ecalDetWindowP_fr_"))
        {
          if (processName == "Scintillation")
          {
            //CreateTree::Instance()->h_photon_2D_receive_Scin->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_r_collect_Scin->Fill(photWL);
            CreateTree::Instance()->SDdetected_rr_S++;
          }
          if (processName == "Cerenkov")
          {
            //CreateTree::Instance()->h_photon_2D_receive_Ceren->Fill(theTrack->GetVertexPosition()[1], theTrack->GetVertexPosition()[2]);
            CreateTree::Instance()->h_phot_lambda_ECAL_r_collect_Ceren->Fill(photWL);
            CreateTree::Instance()->SDdetected_rr_C++;
          }
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }


        if (thePrePVName.contains("world"))
        {
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        else if (thePrePVName.contains("ecalGapP"))
        {
          //std::cout << "in Gap " << thePrePVName << std::endl;
        }
        else if (thePrePVName.contains("wrap"))
        {
          //std::cout << "Crystal " << thePrePVName << std::endl;
        }
        else
        {
          //std::cout << "weird PrePVName " << thePrePVName << std::endl;
        }
        if (thePostPVName.contains("world"))
        {
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        //std::cout << "out of world " << thePrePVName << std::endl;
        G4double tracklength = theStep->GetTrack()->GetTrackLength();
        if (tracklength > maxtracklength)
        {
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
          std::cout << "maximum " << thePrePVName << std::endl;
        }
        if (!propagateCerenkov && (processName == "Cerenkov"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);

        if (!propagateScintillation && (processName == "Scintillation"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);

        if ((nStep == 1) && thePrePVName.contains("ecal"))
        {
          //theTrack->GetDefinition()->GetParticleName()
          //cout<<"position x:"<<global_x<<"  z: "<<global_z<<endl;
          //if (processName == "Scintillation")
          //  CreateTree::Instance()->h_photon_2D_produce_Scin->Fill(global_y, global_z);
          //if (processName == "Cerenkov")
          //  CreateTree::Instance()->h_photon_2D_produce_Ceren->Fill(global_y, global_z);
          if (thePrePVName.contains("ecalCrystalP_f"))
          {
            if (processName == "Scintillation")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_f_produce_Scin->Fill(photWL);
              CreateTree::Instance()->ECAL_f_total_S += 1;
            }
            if (processName == "Cerenkov")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_f_produce_Ceren->Fill(photWL);
              CreateTree::Instance()->ECAL_f_total_C += 1;
            }
          }
          if (thePrePVName.contains("ecalCrystalP_r"))
          {
            CreateTree::Instance()->tot_phot_cer_ECAL_cheren_r_total += 1;
            if (processName == "Scintillation")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_r_produce_Scin->Fill(photWL);
              CreateTree::Instance()->ECAL_r_total_S += 1;
            }
            if (processName == "Cerenkov")
            {
              CreateTree::Instance()->h_phot_lambda_ECAL_r_produce_Ceren->Fill(photWL);
              CreateTree::Instance()->ECAL_r_total_C += 1;
            }
          }
          TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());
          G4int aapdgid = theTrackInfo->GetParentPDGid();

        }

      }
      else
      {
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }
  else
  {

    //count tracks before SCEPCAL at the tracker layers
    if (thePrePVName.contains("world") && thePostPVName.contains("trackerPV_Layer")) // interface between T1 and T2
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePostPVName == Form("trackerPV_Layer%d", iLayer)) // interface between T1 and T2
          CreateTree::Instance()->nTracksTRK[iLayer]++;         //counting tracks crossing the boundary
      }
    }

    //save primary particle position and momentum before timing layer T1 and before ECAL E1
    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_front")) // interface between world and T1
    {

      CreateTree::Instance()->nTracksT1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosT1->at(0) = global_x;
        CreateTree::Instance()->primaryPosT1->at(1) = global_y;
        CreateTree::Instance()->primaryPosT1->at(2) = global_z;

        CreateTree::Instance()->primaryMomT1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomT1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomT1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomT1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if (thePrePVName.contains("world") && thePostPVName.contains("corePV_rear")) // interface between T1 and T2
    {
      CreateTree::Instance()->nTracksT2++; //counting tracks crossing the boundary
    }

    else if ((thePrePVName.contains("world") || thePrePVName.contains("ecalGapP_f") || thePrePVName.contains("ecalDetP_f")) && thePostPVName.contains("ecalCrystalP_f") // interface between world and E1
    )
    {
      CreateTree::Instance()->nTracksE1++; //counting tracks crossing the boundary

      if (theTrack->GetParentID() == 0) // select only the primary particle
      {
        CreateTree::Instance()->primaryPosE1->at(0) = global_x;
        CreateTree::Instance()->primaryPosE1->at(1) = global_y;
        CreateTree::Instance()->primaryPosE1->at(2) = global_z;

        CreateTree::Instance()->primaryMomE1->at(0) = thePrePoint->GetMomentum().x() / GeV;
        CreateTree::Instance()->primaryMomE1->at(1) = thePrePoint->GetMomentum().y() / GeV;
        CreateTree::Instance()->primaryMomE1->at(2) = thePrePoint->GetMomentum().z() / GeV;
        CreateTree::Instance()->primaryMomE1->at(3) = thePrePoint->GetTotalEnergy() / GeV;
      }
    }

    else if ((thePrePVName.contains("ecalCrystalP_f") || thePrePVName.contains("world")) && thePostPVName.contains("ecalCrystalP_r")) // interface between E1 and E2
    {
      CreateTree::Instance()->nTracksE2++; //counting tracks crossing the boundary
    }

    //tracker
    if (thePrePVName.contains("trackerPV_Layer"))
    {
      for (int iLayer = 0; iLayer < 6; iLayer++)
      {
        if (thePrePVName == Form("trackerPV_Layer%d", iLayer))
          CreateTree::Instance()->Edep_Tracker_layer[iLayer] += energy / GeV;
      }
    }

    //timing
    if (thePrePVName.contains("corePV_front"))
    {
      CreateTree::Instance()->depositedEnergyTiming_f += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_f += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_f += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_front_%d", iBar))
          CreateTree::Instance()->Edep_Timing_f_ch[iBar] += energy / GeV;
      }
    }
    if (thePrePVName.contains("corePV_rear"))
    {
      CreateTree::Instance()->depositedEnergyTiming_r += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTiming_r += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTiming_r += energyElec / GeV;
      for (int iBar = 0; iBar < 18; iBar++)
      {
        if (thePrePVName == Form("corePV_rear_%d", iBar))
          CreateTree::Instance()->Edep_Timing_r_ch[iBar] += energy / GeV;
      }
    }

    //ecal

    //hcal
    if (thePrePVName.contains("hcalTile_layer"))
    {
      CreateTree::Instance()->depositedEnergyHCALAct += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALAct += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALAct += energyElec / GeV;
    }
    if (thePrePVName.contains("hcalAbs"))
    {
      CreateTree::Instance()->depositedEnergyHCALPas += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyHCALPas += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyHCALPas += energyElec / GeV;
    }

    if (thePrePVName.contains("world"))
    {
      CreateTree::Instance()->depositedEnergyWorld += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyWorld += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyWorld += energyElec / GeV;
    }

    if (thePrePVName.contains("services"))
    {
      CreateTree::Instance()->depositedEnergyServices += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyServices += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyServices += energyElec / GeV;
    }

    if (thePrePVName.contains("TimingGap"))
    {
      CreateTree::Instance()->depositedEnergyTimingGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyTimingGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyTimingGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalGap"))
    {
      CreateTree::Instance()->depositedEnergyEcalGap += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalGap += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalGap += energyElec / GeV;
    }

    if (thePrePVName.contains("ecalDet"))
    {
      CreateTree::Instance()->depositedEnergyEcalDet += energy / GeV;
      CreateTree::Instance()->depositedIonEnergyEcalDet += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergyEcalDet += energyElec / GeV;
    }

    if (thePrePVName.contains("solenoid"))
    {
      CreateTree::Instance()->depositedEnergySolenoid += energy / GeV;
      CreateTree::Instance()->depositedIonEnergySolenoid += energyIon / GeV;
      CreateTree::Instance()->depositedElecEnergySolenoid += energyElec / GeV;
    }

    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon

  return;
}

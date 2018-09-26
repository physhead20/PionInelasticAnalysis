#define InElastic_v2_cxx
#include "InElastic_v2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <iostream>


// =================================== Other Final State Particle Histograms =================================
TH1D *hNPion = new TH1D("hNPion", "N Pions in InElastic Interaction", 100, 0, 50);
TH1D *hPionDaughterInitalKE = new TH1D("hPionDaughterInitalKE", "Pion Daughter Initial Kinetic Energy", 500, 0, 1000);

TH1D *hNNeutron = new TH1D("hNNeutron", "N Neutrons in InElastic Interaction", 100, 0, 50);
TH1D *hNeutronDaughterInitalKE = new TH1D("hNeutronDaughterInitalKE", "Neutron Daughter Initial Kinetic Energy", 500, 0, 1000);

TH1D *hNPhotons = new TH1D("hNPhotons", "N Photons in InElastic Interaction", 100, 0, 50);
TH1D *hPhotonDaughterInitalKE = new TH1D("hPhotonDaughterInitalKE", "Photon Daughter Initial Kinetic Energy", 500, 0, 1000);

TH1D *hNAtoms = new TH1D("hNAtoms", "N Atoms in InElastic Interaction", 100, 0, 50);
TH1D *hAtomDaughterInitalKE = new TH1D("hAtomDaughterInitalKE", "Atom Daughter Initial Kinetic Energy", 500, 0, 1000);



// ==========================================   Proton Histograms    =========================================
TH1D *hNProton = new TH1D("hNProton", "N Protons in InElastic Interaction", 100, 0, 50);
TH1D *hProtonInitalKE = new TH1D("hProtonInitalKE", "Proton Initial Kinetic Energy", 500, 0, 1000);
TH1D *hNProtonwEngThreshold = new TH1D("hNProtonwEngThreshold", "N Protons w/ Energy threshold", 100, 0, 50);




// ===========================================   Pion Histograms   ===========================================
TH1D *hNPionDaughters = new TH1D("hNPionDaughters", "Number of Daughters in InElastic Interaction", 50, 0, 25);
TH1D *hPionInElasticKE = new TH1D("hPionInElasticKE", "InElastic Pion Initial Kinetic Energy", 40, 0, 1000);
TH1D *hPionDaugtherPDG = new TH1D("hPionDaughterPDG", "Daughter of InElastic Interactions", 20, 0, 10);


// ========================================== Reconstructed Quantities ========================================
TH1D *hRecoNTrack = new TH1D("hRecoNTrack", "Number of Reconstructed Tracks", 50, 0, 25);
TH1D *hRecoDistBewtnMtchTrkandOtherStartPt = new TH1D("hRecoDistBewtnMtchTrkandOtherStartPt", "Distance between end of WC2TPC Trk and Start point of other tracks", 200, 0, 100);
TH1D *hRecoDistBewtnMtchTrkandOtherEndPt = new TH1D("hRecoDistBewtnMtchTrkandOtherEndPt", "Distance between end of WC2TPC Trk and End point of other tracks", 200, 0, 100);


void InElastic_v2::Loop()
{
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;

// ########################
// ### Global Variables ###
// ########################

// --- Counters ---
int ntotalEvents = 0, nPionInElastic = 0, nInFiducialVolume = 0;
int nEnoughChargedPartAtInteraction = 0, nRecoWC2TPCEvents = 0;

float proton_mass = 938.3;
float pion_mass = 139.57; //<---Mass of Pion in MeV
float neutron_mass = 939.5;
float photon_mass = 0.;
float atom_mass = 18629.8;

// ########################################################################
// ### Fiducial Boundry Cuts (used to determine if a track is stopping) ###
// ########################################################################
double XLowerFid = 8;
double XUpperFid = 39;

double YLowerFid = -12;
double YUpperFid = 12;

double ZLowerFid = 8;
double ZUpperFid = 82;




// ################################################
// ### Truth Threshold for proton consideration ###
// ################################################
float ParticleEnergyThreshold = 50;//<---(units of MeV)

// ##############################################
// ###               Event Loop               ###
// ##############################################
//for (Long64_t jentry=0; jentry<nentries;jentry++)
for (Long64_t jentry=0; jentry<15000;jentry++)
   {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   
   
   // #################################
   // ### Proton Momentum variables ###
   // ################################# 
   float g4Proton_Px[100] = {0.}, g4Proton_Py[100] = {0.}, g4Proton_Pz[100] = {0.};
   float g4Proton_momentum[100] = {0.}, g4Proton_kineticEnergy[100] = {0.};
   
   ntotalEvents++;
   
   // === Outputting every nEvents to the screen ===
   if(ntotalEvents % 100 == 0){std::cout<<"Event = "<<ntotalEvents<<std::endl;}
   
   
   
//------------------------------------------------------------------------------------------------------------------   
   // ##################################################################
   // ### First identify if this event is a pion inelastic collision ###
   // ##################################################################
   
   // ### Flag for Pion InElastic Interaction ###
   bool PionInElastic = false;
   
   std::string pi = "pi-Inelastic";
   std::string primary = "primary";
   
   // ### Counter for the iFinal loop ###
   int counter_iFinal = 0;
   
   // ### Variables for the inElastic interaction point ###
   float Interaction_X = -999., Interaction_Y = -999., Interaction_Z = -999.;
   float pionPx = -999., pionPy = -999., pionPz = -999., pionMomentum = -999, PionKineticEnergy = -999;
   float pionTrackID = -999., pionNDaughters = -999.;
   
   // ################################################################
   // ###  Loop over all the particles and identify if the primary ###
   // ### is a pion which ends its life as a inelastic interaction ###
   // ################################################################
   for(auto iFinal = G4FinalProcess->begin(); iFinal != G4FinalProcess->end(); iFinal++)
      {
      if(*iFinal == pi && process_primary[counter_iFinal] == 0)
         {
	 PionInElastic = true;
	 Interaction_X = EndPointx[counter_iFinal];
	 Interaction_Y = EndPointy[counter_iFinal];
	 Interaction_Z = EndPointz[counter_iFinal];
	 pionPx = Px[counter_iFinal] * 1000; //<---Converting to MeV
	 pionPy = Py[counter_iFinal] * 1000; //<---Converting to MeV
	 pionPz = Pz[counter_iFinal] * 1000; //<---Converting to MeV
	 pionMomentum = sqrt( (pionPx*pionPx) + (pionPy*pionPy) + (pionPz*pionPz) );
	 // ###   Calculating the initial Kinetic Energy    ###
         // ### KE = Energy - mass = (p^2 + m^2)^1/2 - mass ###
         PionKineticEnergy = pow( (pionMomentum*pionMomentum) + (pion_mass*pion_mass) ,0.5) - pion_mass;
	 
	 pionTrackID = TrackId[counter_iFinal];
	 pionNDaughters = NumberDaughters[counter_iFinal];
	 }
      
      counter_iFinal++;
      if(PionInElastic){break;}
      }//<---End iFinal loop
   
   // ### Skip the event if it isn't pion in-elastic ###
   if(!PionInElastic){continue;}
   // ### bump the counter ###
   nPionInElastic++;   
//------------------------------------------------------------------------------------------------------------------
   


//------------------------------------------------------------------------------------------------------------------   
   // ##################################################################
   // ### Require that the interaction is within a fiducial boundary ###
   // ##################################################################
   
   // ### Flag for the fiducial boundary ###
   bool InFidVolume = false;
   
   if(Interaction_X > XLowerFid && Interaction_X < XUpperFid &&
      Interaction_Y > YLowerFid && Interaction_Y < YUpperFid &&
      Interaction_Z > ZLowerFid && Interaction_Z < ZUpperFid)
      	{InFidVolume = true;}
   
   // ### Skip the event if it isn't in the fiducial boundary ###
   if(!InFidVolume){continue;}
   // ### bump the counter ###
   nInFiducialVolume++;
    
   // ### Fill the histogram for the InElastic Pions in the fid. volume ###
   hPionInElasticKE->Fill(PionKineticEnergy);
   hNPionDaughters->Fill(pionNDaughters);
//------------------------------------------------------------------------------------------------------------------   


//------------------------------------------------------------------------------------------------------------------
// ####################################################################################
// ### Now we will look at the particles which come out of the in-elastic collision ###
// ####################################################################################  
   
   // ### Setting the flag for enough charged particles at the interaction ###
   bool EnoughChargeParticle = false;
   
   // ### Counter for the iProcess loop ###
   int counter_iProcess = 0;
   
   // ### Proton Variables ###
   int nProtons = 0, nProEThresh = 0;
   
   // ### Other Particle Variables ###
   int nPions = 0, nNeutrons = 0, nPhotons = 0, nAtoms = 0;
   
   int nPiEThresh = 0;
   
   // #######################################
   // ### Looping over the Geant4 process ###
   // #######################################
   for(auto iProcess = G4Process->begin(); iProcess != G4Process->end(); iProcess++)
      {
      
      //std::cout<<" Mother ID = "<<Mother[counter_iProcess]<<", Pion Track ID = "<<pionTrackID<<std::endl;
      // ### Skipping any particle who's mother isn't the primary pion ###
      if(Mother[counter_iProcess] != pionTrackID){counter_iProcess++; continue;}
      //std::cout<<"Process = "<<*iProcess<<", PDG = "<<pdg[counter_iProcess]<<std::endl;
      
      // ### Particles Momentum ###
      float TempPx = Px[counter_iProcess] * 1000; //<---Converting to MeV
      float TempPy = Py[counter_iProcess] * 1000; //<---Converting to MeV
      float TempPz = Pz[counter_iProcess] * 1000; //<---Converting to MeV
      
      // ### Total Momentum ###
      float total_momentum = sqrt( (TempPx*TempPx) + (TempPy*TempPy) + (TempPz*TempPz) );
      
      // ### Classifying the daughters for histogramming ###
      // 0 == Pion
      if(abs(pdg[counter_iProcess]) == 211)
         {
	 hPionDaugtherPDG->Fill(0); 
	 nPions++;
	 
	 float pionKE = pow((total_momentum*total_momentum) + (pion_mass*pion_mass) ,0.5) - pion_mass;
	 hPionDaughterInitalKE->Fill(pionKE);
	 
	 // ###########################################
	 // ### Protons that are above KE threshold ###
	 // ###########################################
	 if(pionKE > ParticleEnergyThreshold)
	    {
	    nPiEThresh++;
	    }//<---end energy threshold
	 }
      // 1 == Proton
      if(abs(pdg[counter_iProcess]) == 2112)
         {
	 hPionDaugtherPDG->Fill(1);
	 }
      // 2 == Neutron
      if(abs(pdg[counter_iProcess]) == 2212)
         {hPionDaugtherPDG->Fill(2); 
	 nNeutrons++;
	 
	 float NeutronKE = pow((total_momentum*total_momentum) + (neutron_mass*neutron_mass) ,0.5) - neutron_mass;
	 hNeutronDaughterInitalKE->Fill(NeutronKE);
	 }
      // 3 == Photon
      if(abs(pdg[counter_iProcess]) == 22)
         {
	 hPionDaugtherPDG->Fill(3); 
	 nPhotons++;
	 
	 float PhotonKE = pow((total_momentum*total_momentum) + (photon_mass*photon_mass) ,0.5) - photon_mass;
	 hPhotonDaughterInitalKE->Fill(PhotonKE);
	 }
      // 4 == Atom
      if(abs(pdg[counter_iProcess]) > 9999)
         {
	 hPionDaugtherPDG->Fill(4); 
	 
	 float AtomKE = pow((total_momentum*total_momentum) + (atom_mass*atom_mass) ,0.5) - atom_mass;
	 hAtomDaughterInitalKE->Fill(AtomKE);
	 
	 nAtoms++;}
      
      // ### Focusing on the protons coming from the interaction ###
      if(pdg[counter_iProcess] == 2212)
         {
	 
	 // ### Getting the Proton Momentum ###
	 g4Proton_Px[nProtons] = Px[counter_iProcess] * 1000; //<---Converting to MeV
	 g4Proton_Py[nProtons] = Py[counter_iProcess] * 1000; //<---Converting to MeV
	 g4Proton_Pz[nProtons] = Pz[counter_iProcess] * 1000; //<---Converting to MeV

	 // ################################################
         // ### Calculating the momentum from the proton ###
         // ################################################
         g4Proton_momentum[nProtons] = sqrt( (g4Proton_Px[nProtons]*g4Proton_Px[nProtons]) + 
	                        (g4Proton_Py[nProtons]*g4Proton_Py[nProtons]) + 
				(g4Proton_Pz[nProtons]*g4Proton_Pz[nProtons]) );
         
   
         // ###   Calculating the initial Kinetic Energy    ###
         // ### KE = Energy - mass = (p^2 + m^2)^1/2 - mass ###
         g4Proton_kineticEnergy[nProtons] = pow( (g4Proton_momentum[nProtons]*g4Proton_momentum[nProtons]) + 
	                                         (proton_mass*proton_mass) ,0.5) - proton_mass;
	 
	 hProtonInitalKE->Fill(g4Proton_kineticEnergy[nProtons]);
	 // ###########################################
	 // ### Protons that are above KE threshold ###
	 // ###########################################
	 if(g4Proton_kineticEnergy[nProtons] > ParticleEnergyThreshold)
	    {
	    nProEThresh++;
	    }//<---end energy threshold
	 
	 
	 nProtons++;
	 }//<---End looking at protons
      
      
      counter_iProcess++;
      }//<---end iProcess Loop
   
   // ### Filling histograms for the final state particles ###
   hNProton->Fill(nProtons);
   hNProtonwEngThreshold->Fill(nProEThresh); 
   
   hNPion->Fill(nPions);  
   hNNeutron->Fill(nNeutrons);
   hNPhotons->Fill(nPhotons);
   hNAtoms->Fill(nAtoms);

   // #########################################################################
   // ### Checking to see if there are enough charged final state daughters ###
   // ###    at the vertex above the energy threshold to tag this event     ###
   // #########################################################################
   if(nProEThresh > 1 || (nProEThresh > 0 && nPiEThresh > 0)) {EnoughChargeParticle = true;}
   
   if(!EnoughChargeParticle){continue;}
   nEnoughChargedPartAtInteraction++;
//------------------------------------------------------------------------------------------------------------------   



//------------------------------------------------------------------------------------------------------------------

   // ###########################################
   // ### Looking at reconstructed quantities ###
   // ###########################################
   
   // ### Checking the number of reconstructed TPC Tracks ###
   hRecoNTrack->Fill(ntracks_reco);
   
   // ### Setting a flag for at least 1 WC2TPC Match
   bool MatchWC2TPC = false;
   
   float WCMatchedTrack_EndX = -999, WCMatchedTrack_EndY = -999, WCMatchedTrack_EndZ = -999;
   
   // ######################################
   // ### Loop over reconstructed tracks ###
   // ######################################
   for(int itrk = 0; itrk < ntracks_reco; itrk++)
      {
      // ### Looking for match = 0 ###
      if(trkWCtoTPCMatch[itrk] == 0)
         {
	 MatchWC2TPC = true;
	 
	 // ### Saving the matched track end points ###
	 WCMatchedTrack_EndX = trkendx[itrk];
	 WCMatchedTrack_EndY = trkendy[itrk];
	 WCMatchedTrack_EndZ = trkendz[itrk];
	 
	 
	 }//<---End WC2TPC Matched Track

      
      }//<---end itrk loop

   // ###############################################################
   // ### Skip any event that doesn't have a matched WC2TPC track ###
   // ###############################################################
   if(!MatchWC2TPC){continue;}
   nRecoWC2TPCEvents++;
//------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------

   // ### Calculating the distance between the matched  ###
   // ### WC2TPC track and all other tracks start point ###
   
   // ### Loop over reconstructed tracks ###
   for(int itrk = 0; itrk < ntracks_reco; itrk++)
      {
      
      // ### Here, we want to skip the matched track ###
      if(trkWCtoTPCMatch[itrk] == 0){continue;}
      
      float deltaRStart = sqrt( pow(WCMatchedTrack_EndX - trkvtxx[itrk], 2 ) + 
                                pow(WCMatchedTrack_EndY - trkvtxy[itrk], 2 ) +
		                pow(WCMatchedTrack_EndZ - trkvtxz[itrk], 2 ));
				
				
      float deltaREnd   = sqrt( pow(WCMatchedTrack_EndX - trkendx[itrk], 2 ) + 
                                pow(WCMatchedTrack_EndY - trkendy[itrk], 2 ) +
		                pow(WCMatchedTrack_EndZ - trkendz[itrk], 2 ));
      
      hRecoDistBewtnMtchTrkandOtherStartPt->Fill(deltaRStart);
      hRecoDistBewtnMtchTrkandOtherEndPt->Fill(deltaREnd);
      }//<---end itrk loop
      
  
   }//<---End jentry loop


// ---------------------------------------------------------------------------------------------------------
std::cout<<std::endl;
std::cout<<"Number of Events 		= "				<<		ntotalEvents		<<std::endl;
std::cout<<"Number of InElastic Events	= "				<<		nPionInElastic		<<std::endl;
std::cout<<"Number interactions in Fid. Vol 	= "			<<	nInFiducialVolume		<<std::endl; 
std::cout<<"Number with >=2 charged particles above energy threshold = "<<	nEnoughChargedPartAtInteraction	<<std::endl;
std::cout<<"Number of events w/ WC2TPC Match   = "			<<	nRecoWC2TPCEvents		<<std::endl;
std::cout<<std::endl;



// ----------------------------------------------------------------------------------------------------------



// ###############################################
// ### Creating a file to output my histograms ###
// ###############################################

TFile myfile("./ProtonInelasticStudy_histo.root","RECREATE");
hNProton->Write();
hNProtonwEngThreshold->Write();
hProtonInitalKE->Write();
hPionInElasticKE->Write();
hNPionDaughters->Write();
hPionDaugtherPDG->Write();
hNPion->Write();  
hNNeutron->Write();
hNPhotons->Write();
hNAtoms->Write();
hPionDaughterInitalKE->Write();
hNeutronDaughterInitalKE->Write();
hPhotonDaughterInitalKE->Write();
hAtomDaughterInitalKE->Write();
hRecoNTrack->Write();
hRecoDistBewtnMtchTrkandOtherStartPt->Write();
hRecoDistBewtnMtchTrkandOtherEndPt->Write();


}//<---End File

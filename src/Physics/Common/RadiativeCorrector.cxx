//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2007 - CA
   Added code to generate the interaction according to a realistic nuclear
   density and made that the default setting.
 @ Dec 01, 2007 - CA
   For COH and ve- interactions setting the vertex on the nuclear boundary
   rather than inside the nucleus.
 @ Sep 15, 2009 - CA
   IsFake() and IsNucleus() are no longer available in GHepParticle. 
   Use pdg::IsPseudoParticle() and pdg::IsIon().
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.
 @ Jun 03, 2016 - JJ (SD)
   Move code to generate the position to a public method, so this code
   can easily be reused by other classes. (Specifically LwlynSmithQELCCPXSec
   and NievesQELCCPXSec to generate a position before calculating the xsec
   when making splines).
*/
//____________________________________________________________________________
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <TMath.h>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TParticlePDG.h>
#include <TF1.h>
#include <TCanvas.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Common/RadiativeCorrector.h"
#include "Physics/Decay/UnstableParticleDecayer.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/RunOpt.h"

using std::count;
using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::controls;
using namespace genie::constants;

//___________________________________________________________________________
RadiativeCorrector::RadiativeCorrector() :
EventRecordVisitorI("genie::RadiativeCorrector")
{
  fInitState = 0;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
RadiativeCorrector::RadiativeCorrector(string config) :
EventRecordVisitorI("genie::RadiativeCorrector", config)
{
}
//___________________________________________________________________________
RadiativeCorrector::~RadiativeCorrector()
{
  if (fInitState)        delete fInitState;
  //rad_gOptRunNu = 0;
  //rad_kDefOptNtpFormat = kNFGHEP;
}
//___________________________________________________________________________
void RadiativeCorrector::BuildInitialState(const InitialState & init_state)
{

  if(fInitState) delete fInitState;
  fInitState = new InitialState(init_state);

  //this->AssertIsValidInitState();
}
//___________________________________________________________________________
void RadiativeCorrector::ProcessEventRecord(GHepRecord * evrec) const 
{
  LOG("RadiativeCorrector", pDEBUG) << "Welcome to the Radiative corrector using the "<<fModel<<" model";
  if (fISR) LOG("RadiativeCorrector", pDEBUG) <<"This is ISR";
  else LOG("RadiativeCorrector", pDEBUG) <<"This is FSR";
  // decay a photon with dE from the electron
  evrec->SetPrintLevel(13); 
  //evrec->Print();
  //std::cout<<"\n";
  bool radDone = false;
 
  Interaction * interaction = evrec->Summary();
  Kinematics * kine = interaction->KinePtr();
  InitialState * init_state_ptr = interaction -> InitStatePtr();
  Target * target_ptr = init_state_ptr -> TgtPtr();
  int Z = target_ptr->Z();
  
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  unsigned int ipos = 0;
  //unsigned int hitNucleonPosition; 
  while( (p = (GHepParticle *) piter.Next()) )
  {
     if( this->ToBeDecayed(p)) {
        TLorentzVector p4 = *(p->P4());
        TLorentzVector x4 = *(p->X4());
	double e = p4.E();
	double energyLoss = 0.;
        double e_gamma_max;
	double e_gamma_min = 1E-25;
        
        e_gamma_max = 0.2*p4.E();
        
	double L1,L2,b,lambda_e,g,power_hi,power_lo;
	if (fModel =="simc" || fModel == "simple") {
 
	   if (Z==1) {
	      L1 = 5.31;
	      L2 = 6.144; 
	   }
	   else {
              L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
	      L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
	   }
           b = (1./9)*(12 + float(Z+1)/(Z*L1 + L2));
	}
        if (fISR) lambda_e = (kAem/kPi)*( 2*TMath::Log(2*p->P4()->P()/kElectronMass) -1 + TMath::Log(0.5*(1-fP4l.CosTheta())) );//+ 2*TMath::Log(init_state_ptr->GetProbeP4(kRfLab)->P()/fP4l.E()) + TMath::Log(0.5*(1-fP4l.CosTheta() ) ) );
	else lambda_e =      (kAem/kPi)*( 2*TMath::Log(2*p->P4()->P()/kElectronMass) -1 + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta())) );//+ 2*TMath::Log(init_state_ptr->GetProbeP4(kRfLab)->P()/kine->FSLeptonP4().E()) + TMath::Log(0.5*(1-kine->FSLeptonP4().CosTheta() ) ) );
        g = b*fThickness + lambda_e;

	power_hi = pow(e_gamma_max,g);
	power_lo  = pow(fCutoff,g);
	TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
        f->SetParameter(0,g);
        f->SetParameter(1,power_hi - power_lo);
        energyLoss = f->GetRandom();
	delete f;
        if (energyLoss<fCutoff) energyLoss = 0.; 

	double momentumLoss = energyLoss;
	double ptLoss;
	double pzLoss;
	if (p4.Pz()==p4.E()) // for the z direction going probe theta = -nan 
	{
	   ptLoss = 0.;
	   pzLoss = momentumLoss;
	}
	else
	{
	   ptLoss = momentumLoss*sin(p4.Theta()); 
	   pzLoss = momentumLoss*cos(p4.Theta());
	}
	TLorentzVector p4RadGamma;
	p4RadGamma.SetPxPyPzE(ptLoss*cos(p4.Phi()),ptLoss*sin(p4.Phi()),pzLoss,energyLoss);
	TLorentzVector p4tag = p4 - p4RadGamma;

	if (fISR && !radDone && energyLoss>0) {
	  init_state_ptr->SetProbeP4(p4tag);
	  LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG=" << p->Pdg();
          evrec->AddParticle(p->Pdg(), kIStCorrectedProbe, ipos,-1,-1,-1, p4tag, x4);
          LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG= 22";
          evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
	  radDone = true;
	}
	
        if (!fISR && !radDone) {
 	   //-- Update the event weight for each weighted particle decay
	   float radcor_weight = 1.;
	   if (energyLoss>0) {
	       LOG("RadiativeCorrector", pDEBUG) << "performing FSR correction for: " << p->Name();
               //-- Mark it as a 'decayed state' & add its daughter links
       	       p->SetStatus(kIStDecayedState);
               ////-- Add the mom & daughters to the event record
               LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG=" << p->Pdg();
	       evrec->AddParticle(p->Pdg(), kIStStableFinalState, ipos,-1,-1,-1, p4tag, x4);
               LOG("RadiativeCorrector", pDEBUG) << "Adding daughter... PDG= 22";
	       evrec->AddParticle(22, kIStStableFinalState, ipos,-1,-1,-1, p4RadGamma, x4);
               ////-- Printouts
	       //evrec->Print();
	       //std::cout<<"\n";
	       radDone = true;
	   }
	}
     }
     //if (fISR && (p->Status() ==1)) p->SetStatus(kIStIntermediateState);
     if (fISR && (p->Status() ==1) && (p->Pdg()==11)) {
        
	p->SetStatus(kIStDecayedState);
	LOG("RadiativeCorrector", pDEBUG) <<"After setting the particles status to decyed and momentum of the nucleon to be zero, pdg "<<p->Pdg()<<" mass "<<p->Mass();
        //evrec->Print();
	//std::cout<<"\n";
        //kine->ClearRunningValues();
        //kine->Reset();
     }
     ipos++;
  } // loop over particles
  //LOG("RadiativeCorrector", pINFO) << "The final state lepton energy is "<<evrec->FinalStatePrimaryLepton()->P4()->E();


}
//___________________________________________________________________________
bool RadiativeCorrector::ToBeDecayed(GHepParticle * particle) const
{
   LOG("RadiativeCorrector", pDEBUG) <<"Particle ToBeDecayed pdg "<<particle->Pdg()<<" status "<<particle->Status()<<" mother "<<particle->FirstMother();
   if(particle->Pdg() != 0) {
     if (fISR) {
     	if ((particle->Status() == 0) && (pdg::IsElectron(particle->Pdg()) || pdg::IsPositron(particle->Pdg()))) return true;
     	else return false;	
     }
     else {
     	if ((particle->Status() == 1) && (pdg::IsElectron(particle->Pdg()) || pdg::IsPositron(particle->Pdg()))) return true;
     	else return false;	
     }
   }
   return false;
}
//___________________________________________________________________________
bool RadiativeCorrector::ParticleWasRadiated(GHepParticle * particle) const
{
  std::cout<<"particle was radiated pdg "<<particle->Pdg()<<" status "<<particle->Status()<<" mother "<<particle->FirstMother() <<std::endl;
  if ( (particle->Pdg() ==11) && ((particle->Status() == 5) || (particle->Status() ==3))) return true;
  else return false;
}
//___________________________________________________________________________
void RadiativeCorrector::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeCorrector::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RadiativeCorrector::LoadConfig(void)
{
  GetParam( "ISR",fISR);
  GetParam( "RadiativeCorrectionModel" , fModel);
  GetParam( "RadiativeCorrectionThickness_1000010010",fThickness_1000010010);
  GetParam( "RadiativeCorrectionThickness_1000020040",fThickness_1000020040);
  GetParam( "RadiativeCorrectionThickness_1000060120",fThickness_1000060120);
  GetParam( "RadiativeCorrectionThickness_1000260560",fThickness_1000260560);
  GetParam( "RadiativeCorrectionCutoff",fCutoff);
  GetParam( "RadiativeCorrectionDoInternal",fDoInternal);
  
  /*
  const std::string keyStart = "Thickness@Pdg=";

  RgIMap entries = GetConfig().GetItemMap();

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it){
    const std::string& key = it->first;
    int pdg = 0;
    int Z = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str())) {
      pdg = atoi(key.c_str() + keyStart.size());
      Z = pdg::IonPdgCodeToZ(pdg);
    }
    if (0 != pdg && 0 != Z) {
      ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double thickness ;
      GetParam( rgkey, thickness ) ;
      thickness = TMath::Max(thickness, 0.);
      LOG("RadiativeCorrector", pINFO) << "Nucleus: " << pdg << " -> using Thickness =  " << thickness << " radiation lengths";
      fThicknesses.insert(map<int,double>::value_type(Z,thickness));
    }
  }*/
}
//____________________________________________________________________________
void RadiativeCorrector::SetISR(bool isr)
{
  fISR = isr;
}
//____________________________________________________________________________
void RadiativeCorrector::SetModel(std::string model)
{
  fModel = model;
}
//____________________________________________________________________________
void RadiativeCorrector::SetQ2(double Q2)
{
  fQ2 = Q2;
}
//____________________________________________________________________________
void RadiativeCorrector::SetP4l(TLorentzVector p4l)
{
  fP4l = p4l;
}
//____________________________________________________________________________
void RadiativeCorrector::SetCutoff(double cutoff)
{
  fCutoff = cutoff;
}
//____________________________________________________________________________
void RadiativeCorrector::SetThickness(int tgtpdg, double thickness_1000010010, double thickness_1000020040, double thickness_1000060120, double thickness_1000260560)
{
  if      (tgtpdg == 1000010010) fThickness = thickness_1000010010;
  else if (tgtpdg == 1000020040) fThickness = thickness_1000020040;
  else if (tgtpdg == 1000060120) fThickness = thickness_1000060120;
  else if (tgtpdg == 1000260560) fThickness = thickness_1000260560;
}
//____________________________________________________________________________
void RadiativeCorrector::SetDoInternalRad(bool doInternal)
{
  fDoInternal = doInternal;
}
//____________________________________________________________________________
void RadiativeCorrector::Configure(const InitialState & is)
{
  /*InitialState init_state(is.TgtPdg(), is.CorrectProbePdg()); // filter any other init state info

  ostringstream mesg;
  mesg << "Configuring for initial state: `"
       << init_state.AsString();

  LOG("RadiativeCorrector", pNOTICE)
        << utils::print::PrintFramedMesg(mesg.str(), 0, '*');

  this -> BuildInitialState            (init_state);
  //this -> BuildGeneratorList           ();
  //this -> BuildInteractionGeneratorMap ();
  //this -> BuildInteractionSelector     ();

  LOG("RadiativeCorrector", pINFO) << "Done configuring. \n";*/
}
//____________________________________________________________________________
/*void RadiativeCorrector::SetE(double e) 
{
	r_e = e;
}*/
//____________________________________________________________________________

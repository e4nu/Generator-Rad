//____________________________________________________________________________
/*!

\class    genie::RESKinematicsGenerator

\brief    Generates resonance event (v+N->l+Resonance) kinematics.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  November 18, 2004

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _RES_KINEMATICS_GENERATOR_H_
#define _RES_KINEMATICS_GENERATOR_H_

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"

class TF2;

namespace genie {

class RESKinematicsGenerator : public KineGeneratorWithCache {

public :
  RESKinematicsGenerator();
  RESKinematicsGenerator(string config);
 ~RESKinematicsGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig      (void);
  double ComputeMaxXSec  (const Interaction * interaction) const;
  TLorentzVector GetFinalStateLeptonKinematic(GHepRecord * evrec, double E, double gy, double gQ2) const;


  mutable TF2 * fEnvelope; ///< 2-D envelope used for importance sampling
  double fWcut;            ///< Wcut parameter in DIS/RES join scheme
  bool fDoRadiativeCorrection;
  bool fDoInternal;
  std::string fModel;               ///< to distinguish between differnt models, right now simc / vanderhagen
  double      fCutoff;
  double      fThickness;           ///< thicnknesses of targets in CLAS in radiation length 
  double      fThickness_1000010010;           ///< thicnknesses of targets in CLAS in radiation length 
  double      fThickness_1000020040;           ///< thicnknesses of targets in CLAS in radiation length 
  double      fThickness_1000060120;           ///< thicnknesses of targets in CLAS in radiation length 
  double      fThickness_1000260560;           ///< thicnknesses of targets in CLAS in radiation length 
};

}      // genie namespace
#endif // _RES_KINEMATICS_GENERATOR_H_

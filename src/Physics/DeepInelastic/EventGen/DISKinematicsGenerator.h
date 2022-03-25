//____________________________________________________________________________
/*!

\class    genie::DISKinematicsGenerator

\brief    Generates values for the kinematic variables describing DIS v
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Part of its implementation, related with the caching and retrieval of
          previously computed values, is inherited from the KineGeneratorWithCache
          abstract class.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _DIS_KINEMATICS_GENERATOR_H_
#define _DIS_KINEMATICS_GENERATOR_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class DISKinematicsGenerator : public KineGeneratorWithCache {

public :
  DISKinematicsGenerator();
  DISKinematicsGenerator(string config);
  ~DISKinematicsGenerator();

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

#endif // _DIS_KINEMATICS_GENERATOR_H_

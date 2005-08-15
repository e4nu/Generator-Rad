//____________________________________________________________________________
/*!

\class   genie::GHepRecord

\brief   Generated Event Record: STDHEP-like record and Summary Information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_H_
#define _GHEP_RECORD_H_

#include <TClonesArray.h>
#include "Interaction/Interaction.h"

class TLorentzVector;

namespace genie {

class GHepParticle;

class GHepRecord : public TClonesArray {

public :

  GHepRecord();
  GHepRecord(int size);
  GHepRecord(const GHepRecord & record);
  ~GHepRecord();

  //-- methods to attach / get summary information

  virtual Interaction * GetInteraction    (void) const;
  virtual void          AttachInteraction (Interaction * interaction);

  //-- methods to search the GHEP (STDHEP-like) record

  virtual GHepParticle * GetParticle    (int position) const;
  virtual GHepParticle * FindParticle   (int pdg, int status, int start) const;

  virtual int ParticlePosition(int pdg, int status,  int start) const;
  virtual int ParticlePosition(GHepParticle * particle, int start) const;

  //-- operations on the record

  virtual void Copy                    (const GHepRecord & record);
  virtual void ShiftVertex             (const TLorentzVector & vec4);
  virtual void ResetGHepRecord         (void);
  virtual void CompactifyDaughterLists (void);

  //-- provide a simplified wrapper of the 'new with placement'
  //   TClonesArray object insertion method

  //   ALWAYS use these methods to insert new particles as they check
  //   for the compactness of the daughter lists.
  //   Note that the record might be automatically re-arranged as the
  //   result of your GHepParticle insertion

  virtual void AddParticle (const GHepParticle & p);
  virtual void AddParticle (int pdg, int status,
                     int mom1, int mom2, int dau1, int dau2,
                        const TLorentzVector & p, const TLorentzVector & v);
  virtual void AddParticle (int pdg, int status,
                     int mom1, int mom2, int dau1, int dau2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t);

  //-- methods to switch on/off and ask for event record flags

  virtual void SwitchIsPauliBlocked (bool on_off);
  virtual bool IsPauliBlocked (void) const { return fEventIsPauliBlocked; }
  virtual bool IsForbidden    (void) const;

protected:

  // Summary information for the Initial State, Process Type & Kinematics
  Interaction * fInteraction;

  // Frags for the generated event
  bool fEventIsPauliBlocked;

  // Utility methods
  void InitGHepRecord  (void);

  // Methods used by the daughter list compactifier
  virtual void UpdateDaughterLists    (void);
  virtual bool HasCompactDaughterList (int pos);
  virtual void SwapParticles          (int i, int j);
  virtual void FinalizeDaughterLists  (void);
  virtual int  FirstNonInitStateEntry (void);

private:

ClassDef(GHepRecord, 1)

};

}      // genie namespace

#endif // _GHEP_RECORD_H_

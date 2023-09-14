#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"

class DetectorConstruction;

/// Action initialization class.
///

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization();
    virtual ~ActionInitialization();

    //virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif

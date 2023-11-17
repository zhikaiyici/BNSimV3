#ifndef ListAccumulable_h
#define ListAccumulable_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"
#include "globals.hh"
#include <list>

class ListAccumulable :
    public G4VAccumulable
{
public:
    ListAccumulable(const G4String& name) : G4VAccumulable(name/*, 0*/), fMyList() {};
    virtual ~ListAccumulable() {};

    virtual void Merge(const G4VAccumulable& other);
    virtual void Reset();

private:
    std::list<G4double> fMyList;
};

#endif

/*#include "G4VAccumulable.hh"
#include "globals.hh"
#include <map>
class ProcCounterAccumulable : public G4VAccumulable
{
public:
    ProcCounterAccumulable(const G4String& name)
        : G4VAccumulable(name, 0), fProcCounter() {}
    virtual ~ProcCounterAccumulable() {}

    void CountProcesses(G4String procName);

    virtual void Merge(const G4VAccumulable& other);
    virtual void Reset();

private:
    std::map<G4String, G4int> fProcCounter;
};

void ProcCounterAccumulable::Merge(const G4VAccumulable& other)
{
    const ProcCounterAccumulable& otherProcCounterAccumulable
        = static_cast<const ProcCounterAccumulable&>(other);

    std::map<G4String, G4int>::const_iterator it;
    for (it = otherProcCounterAccumulable.fProcCounter.begin();
        it != otherProcCounterAccumulable.fProcCounter.end(); ++it) {

        G4String procName = it->first;
        G4int otherCount = it->second;
        if (fProcCounter.find(procName) == fProcCounter.end()) {
            fProcCounter[procName] = otherCount;
        }
        else {
            fProcCounter[procName] += otherCount;
        }
    }
}

void ProcCounterAccumulable::Reset()
{
    fProcCounter.clear();
}*/

#include "ListAccumulable.hh"
#include <list>

void ListAccumulable::Merge(const G4VAccumulable& other)
{
	const ListAccumulable& otherListAccumulable = static_cast<const ListAccumulable&>(other);
    std::list<G4double> otherList = otherListAccumulable.fMyList;
    fMyList.merge(otherList);
}

void ListAccumulable::Reset()
{
	fMyList.clear();
}

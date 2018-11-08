#include "PGOneI_PerProtein.h"
#include <string>
#include <vector>

PerProtein::PerProtein()
{
  RealProteinIndex = 0;
  LengthOfProtein = 0;
  CurrentPrefixSize = 0;
}


PerProtein::PerProtein(const PerProtein & rval)
{
  RealProteinIndex = rval.RealProteinIndex;
  LengthOfProtein = rval.LengthOfProtein;
  CurrentPrefix = rval.CurrentPrefix;
  CurrentPrefixSize = rval.CurrentPrefixSize;

}

void PerProtein::InitiatePerProtein(const string & OneProtein,const char & FirstPrefix, const int & ProteinIndexInput) {
  RealProteinIndex = ProteinIndexInput;
  LengthOfProtein = OneProtein.size();
  for (int i = 0; i < LengthOfProtein-1; i++ )
    if ( OneProtein[i] == FirstPrefix ) {
      CurrentPrefix.push_back(i);
      CurrentPrefixSize++;
    }
}

void PerProtein::UpdateCurrentPrefix(const string & OneProtein,
                                     const char & NewPrefix) {
    for (int PositionIndex = CurrentPrefixSize-1; PositionIndex >= 0;
         PositionIndex--)
      if (OneProtein[++CurrentPrefix[PositionIndex]] != NewPrefix)
	 DeleteOnePosition(PositionIndex);
  }

void PerProtein::RightMovePointer() {
  for (int PositionIndex = CurrentPrefixSize - 1; PositionIndex >= 0;
       PositionIndex-- )
    if (CurrentPrefix[PositionIndex] >= LengthOfProtein-2)  // think of extra #
      DeleteOnePosition(PositionIndex);  // Current pointer point to end of the sequence;
    else CurrentPrefix[PositionIndex]++;
}

void PerProtein::DeleteOnePosition(const int & PositionIndex) {
  CurrentPrefix[PositionIndex] = CurrentPrefix[--CurrentPrefixSize];
  CurrentPrefix.pop_back();
}

void PerProtein::Copy(const PerProtein & rval) {
  RealProteinIndex = rval.RealProteinIndex;
  LengthOfProtein = rval.LengthOfProtein;
  //EndOfPrefix = rval.EndOfPrefix;
  CurrentPrefix = rval.CurrentPrefix;
  CurrentPrefixSize = rval.CurrentPrefixSize;
}

  //int PerProtein::GetRealProteinIndex() {return RealProteinIndex;}
  int PerProtein::GetNumberOfPositions() {return CurrentPrefixSize;}
  int PerProtein::GetLengthOfProtein() {return LengthOfProtein;}

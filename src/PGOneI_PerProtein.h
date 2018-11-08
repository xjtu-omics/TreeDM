// PGOneI_PerProtein.h

#ifndef PERPROTEIN_H
#define PERPROTEIN_H


#include <string>
#include <vector>

using namespace std;

class PerProtein
{
public:
  PerProtein();
  PerProtein(const PerProtein & rval);

  void InitiatePerProtein(const std::string & OneProtein, const char & FirstPrefix, const int & ProteinIndexInput);
  // For a given protein, find all occurence of one amino acid subtype(FirstPrefix);

  void UpdateCurrentPrefix(const string & OneProtein, const char & NewPrefix);
  // For a given protein, find all occurence of new prefix within the distance of
  // Max_Pat_Length to the start of the motif;
  // If not, delete CurrentPrefix and EndOfPrefix pair;

  void DeleteOnePosition(const int & PositionIndex);
  // When motif cannot grow for current PositionIndex, delete it;

  void Copy(const PerProtein & rval);
  // Make a copy of one PerProtein Obj;

  void RightMovePointer();

  int GetRealProteinIndex();
  int GetNumberOfPositions();
  int GetLengthOfProtein();
  vector<int> CurrentPrefix;
  int RealProteinIndex;


private:
  //int RealProteinIndex;
  int LengthOfProtein;
  //vector<int> EndOfPrefix;
  //vector<int> CurrentPrefix;
  int CurrentPrefixSize;
};

#endif // PERPROTEIN_H

#include "PGOneI_PD.h"
#include "PGOneI_PerProtein.h"

#include <vector>

ProjectedDatabase::ProjectedDatabase()
{
  Support = 0;
  CurrentWildCardLength = 0;
}


ProjectedDatabase::ProjectedDatabase(const ProjectedDatabase & rval,
                                     bool empty) {
  Support = rval.Support;
  CurrentWildCardLength = rval.CurrentWildCardLength;
  Prefix = rval.Prefix;
  if ( ! empty ) ForEachProtein = rval.ForEachProtein;
}

void ProjectedDatabase::InitiateProData(const char & FirstPrefix,
                const vector <string> & SequencesDatabase) {
  for (int ProteinIndex = 0; ProteinIndex < (int)SequencesDatabase.size();
       ProteinIndex++) {
    PerProtein TempPerProtein;
    TempPerProtein.InitiatePerProtein(SequencesDatabase[ProteinIndex],
                                      FirstPrefix, ProteinIndex);
    ForEachProtein.push_back(TempPerProtein);
  }
  Support = ForEachProtein.size();
  Prefix.push_back(FirstPrefix);
}

void ProjectedDatabase::UpdateProData(const ProjectedDatabase & TempProData,
      const char & NewPrefix, const vector <string> & SequencesDatabase) {
  CurrentWildCardLength = 0;
  int PerProteinIndex = TempProData.ForEachProtein.size() - 1;
  while (PerProteinIndex >= 0) {
    PerProtein TempProt = TempProData.ForEachProtein[PerProteinIndex];
    TempProt.UpdateCurrentPrefix(
          SequencesDatabase[TempProt.RealProteinIndex], NewPrefix);
    if (TempProt.GetNumberOfPositions() != 0)
      ForEachProtein.push_back(TempProt);
    PerProteinIndex--;
  }
  Prefix.push_back(NewPrefix);
  Support = ForEachProtein.size();
}

void ProjectedDatabase::MoveAll() {
  CurrentWildCardLength++;
  Prefix.push_back('x');
  int PerProteinIndex = ForEachProtein.size() - 1;
  while (PerProteinIndex >= 0) {
    ForEachProtein[PerProteinIndex].RightMovePointer();
    PerProteinIndex--;
  }
}

void ProjectedDatabase::DeleteOneProtein(const int & PerProteinIndex) {
  int LastProteinIndex = ForEachProtein.size() - 1;
  if (PerProteinIndex != LastProteinIndex)
    ForEachProtein[PerProteinIndex].Copy(ForEachProtein[LastProteinIndex]);
  ForEachProtein.pop_back();
}

void ProjectedDatabase::Copy(const ProjectedDatabase &rval) {
  Support = rval.Support;
  Prefix = rval.Prefix;
  ForEachProtein = rval.ForEachProtein;
}

int ProjectedDatabase::GetSupport() {return Support;}

int ProjectedDatabase::GetPrefixSize() {return Prefix.size();}

int ProjectedDatabase::GetCurrentWildCardLength() {
  return CurrentWildCardLength; }

int ProjectedDatabase::GetTotalWildCardLength() {
  int TempTotalWildCardLength = 0;
  for (int PrefixIndex = 0; PrefixIndex < (int)Prefix.size();
       PrefixIndex++) {
    if (Prefix[PrefixIndex]=='x')
      TempTotalWildCardLength++;
  }
  return TempTotalWildCardLength;
}

vector<char> ProjectedDatabase::GetPrefix(){return Prefix;}

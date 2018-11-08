
#ifndef PROJECTEDDATABASE_H
#define PROJECTEDDATABASE_H


#include <string>
#include <vector>

#include "PGOneI_PerProtein.h"


using namespace std;


class ProjectedDatabase {
public:
  ProjectedDatabase();
  ProjectedDatabase(const ProjectedDatabase & rval, bool empty);

  void InitiateProData(const char & FirstPrefix, const vector <string> & SequencesDatabase);
  // For a given input and one bite prefix, find all occurence of the prefix in each proteins;

  void UpdateProData(const ProjectedDatabase & TempProData,
     const char & NewPrefix, const vector <string> & SequencesDatabase);
  // For new prefix, check whether the right side of the recorded location has that prefix.
  // If yes, record new location; if not, delete that record

  void DeleteOneProtein(const int & PerProteinIndex);
  // If there is no record of occurence in one protein, delete that protein from
  // the projected database

  void Copy(const ProjectedDatabase &rval);

  void MoveAll();

  int GetSupport();
  int GetPrefixSize(); // How long is the pattern.
  int GetCurrentWildCardLength();
  int GetTotalWildCardLength();  // For AxxTxF, it is 2 + 1 = 3
  vector<char> GetPrefix();

  vector<PerProtein> ForEachProtein;
  vector<char> Prefix;

private:
  int Support;  // How many proteins support current prefix
  int CurrentWildCardLength; // For AxxTxF, it is 0; for AxxTxFxx, it is 2
  //vector<char> Prefix;
  //  vector<PerProtein> ForEachProtein;
};

#endif // PROJECTEDDATABASE_H

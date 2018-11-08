#include <stdio.h>
#include <cstdlib>
#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <numeric>

//#include "cmdline.h"
#include "PGOneI_PD.h"
#include "PGOneI_PerProtein.h"



using namespace std;
//const string AASubType = "GATC"; // for DNA sequences
const string AASubType = "GASTCVLIMPFYWDENQHKR";

struct RPIPP {        // RealProteinIndexPatternPosition
	int RPI; //RealProteinIndex
	vector<int> PP;  // PatternPosition
};
//create RPIPP to store current pattern position information
//use OUTRPIPP to store current all_position_informations of one_pattern
struct MotifSupportPosition {
	vector <char> Motif;
	int Support;
	vector<RPIPP> OUTRPIPP;
};


void helpinformation(void) {

	std::cerr << "\nUsage: ./TreeDM [options] ...\n\n"
		<< "options:\n"
		<< "  -n   <int>   The minimal number of non_wildcard residues in the pattern to be reported\n"
		<< "  -s   <float>   The minimal proportion of proteins haveing the pattern\n"
		<< "  -m   <int>   The maximal length of wildcard in the patterns\n"
		<< "  -i   <string>   The pathway and filename of sequences in fasta format\n"
		<< "  -p   <string>   The pathway and filename of the PatternSpace file£¨defaults: ./ppfile£©\n"
		<< "  -o   <string>   The pathway of the distance matrix(defaults: ./)"
		<< "  \n"
		<< "  -h/?   help pritn this message\n"
		<< std::endl;
	exit(1);
}

int non_wildcard;
float support;
int max_wildcard;
string inputfile;
string pattern_space = "./ppfile";
string outputDM = "./";

int tree_options(int aargc, char *aargv[]) {
	int i;
	for (i = 1; i < aargc; i++) {
		if (aargv[i][0] != '-') return i;
		switch (aargv[i][1]) {
		case 'n': non_wildcard = atoi(aargv[++i]); break;
		case 's': support = atof(aargv[++i]); break;
		case 'm': max_wildcard = atoi(aargv[++i]); break;
		case 'i': inputfile = aargv[++i]; break;
		case 'p': pattern_space = aargv[++i]; break;
		case 'o': outputDM = aargv[++i]; break;
		break;
		case 'h': helpinformation();
		case '?': helpinformation();
		}
	}
	return i;
}



//struct PatternPosition {
//	int RealProteinIndex;
//	vector<int> PositionInProtein;
//};
//struct PatternInf {
//	string FrequentPattern;
//	int FrequentPatternSupport;
//	vector<PatternPosition> FrequentPatternPosition;
//};
int MaxPatternLength = 0;
vector<string> AllInputFile;

void mining(ProjectedDatabase & TempProData,
	const vector <string>& SequencesDatabase,
	const double & MinSupRatio,
	const int & Min_Pat_Length,
	const int & Max_WildCard_Length,
	int & MaxPL,
	vector <MotifSupportPosition> & TempOutPut,
	vector <MotifSupportPosition> & OutPut);

vector<vector<int> > BinaryPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut,
	const vector<string> &Label, string OutFileName);
vector<vector<int> > IDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut,
	const vector<string> &Label, string OutFileName);
vector<vector<int> > FIDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut,
	const vector<string> &Label, string OutFileName);
vector<vector<int> > WIDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut,
	const vector<string> &Label);


vector<vector<int> > ConsiderWeightDuplicateOverlapPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut,
	const vector<string> &Label, string OutFileName);

void JaccardDistanceMatrix(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName);
void TanimotoDistanceMatrix(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName);
void JensenShannonDivergence(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, const string OutFileName);
void EuclideanDistanceMatrix(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName);
void CosineSimilarityMatirx(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName);

int main(int argc, char *argv[])
{
	if (argc == 1) helpinformation();
	for (int i = 0; i < argc; i++)
	{
		std::cout << argv[i] << ' ';
	}

	int options = tree_options(argc, argv);

	ifstream  inf_seq(inputfile.c_str());
	if (!inf_seq)
	{
		cerr << "Sorry, cannot find the input_file" << endl;
	}

	ofstream outpatterns(pattern_space.c_str());
	if (!outpatterns)
		cout << "Sorry, cannot write to the patterns_file: " << endl;
	
	
	vector < string > SequencesDatabase;
	vector < string > Label;
	string TempName, TempStr;
	while (inf_seq >> TempName >> TempStr)
	{
		if (!TempName.empty())
			Label.push_back(TempName);
		if (!TempStr.empty())
			SequencesDatabase.push_back(TempStr + "#");
	};
//	for (int i = 0; i <= Label.size() - 1; i++)
//	{
//		for (int j = 0; j <= Label[i].size() - 1; j++)
//		{
//			cout << Label[i][j];
//		}
//		cout << endl;
//	}
	cout << "sequencesdatabase size: " << SequencesDatabase.size() << endl;
//	 ##################### end of input sequences database ###########
//
//	 ##################### prefixspan ################################
	ProjectedDatabase EmptyProData;
	vector <MotifSupportPosition> OutPut;
	vector <MotifSupportPosition>  TempOutPut;
	MotifSupportPosition TestMotif;
	vector<char> FirstMotif;
	FirstMotif.push_back('@');
	TestMotif.Motif = FirstMotif;
	TestMotif.Support = SequencesDatabase.size();
	TempOutPut.push_back(TestMotif);
	
	
	//vector <MotifSupportPosition> OutPut;
	cout << "begin of mining run!" << endl;
	for (int AAIndex = 0; AAIndex < 20; AAIndex++)
	{
		ProjectedDatabase ProData(EmptyProData, false);// create projecteddatabase
		ProData.InitiateProData(AASubType[AAIndex], SequencesDatabase);
		cout << "AA is ..." << AASubType[AAIndex] << endl;
		mining(ProData, SequencesDatabase, support, non_wildcard, max_wildcard, MaxPatternLength, TempOutPut, OutPut);
	}
	if (TempOutPut[TempOutPut.size() - 1].Motif[0] != '@')
	{
		OutPut.push_back(TempOutPut[TempOutPut.size() - 1]);// store the last one
	}
	TempOutPut[0] = TempOutPut[TempOutPut.size() - 1];// delete first element
	TempOutPut.pop_back();
	
	cout << "the maxpl is " << MaxPatternLength << endl;


	cout << "End of mining run!" << endl;
	cout << "close pattern_numbers = " << OutPut.size() << endl;
	cout << "all pattern numbers is " << TempOutPut.size() << endl;

 	outpatterns << "The numbers of FP is " << OutPut.size() << endl;
	for (int OutPutIndex = 0; OutPutIndex < (int)OutPut.size(); OutPutIndex++)
	{
		for (int PrefixIndex = 0; PrefixIndex < (int)OutPut[OutPutIndex].Motif.size(); PrefixIndex++)
		{
			outpatterns << OutPut[OutPutIndex].Motif[PrefixIndex];
		}
		outpatterns << "\t" << OutPut[OutPutIndex].Support << endl/* "\t"*/;
		//kyy*******************************************************pattern position******************************************************
		for (int preProindex = 0; preProindex <= OutPut[OutPutIndex].OUTRPIPP.size() - 1; preProindex++)
		{
			outpatterns << "(" << OutPut[OutPutIndex].OUTRPIPP[preProindex].RPI << "|";
			for (int positionindex = 0; positionindex <= OutPut[OutPutIndex].OUTRPIPP[preProindex].PP.size() - 1; positionindex++)
			{
				outpatterns << OutPut[OutPutIndex].OUTRPIPP[preProindex].PP[positionindex] << "\t";
			}
			outpatterns << ")";
		}
		outpatterns << "\n";
		//	//kyy******************************************************pattern position********************************************************
	}
	//************OutPut***********
    //string DMPath;
    //DMPath.assign(argv[6],strlen(argv[6]));

	if (OutPut.size() == 0)
	{
			cout << "no distance matrix file" << endl;
			return 0;
	}
	else
	{

		vector<vector<int> > TotalPatternVector;





		TotalPatternVector = WIDFWeightPatternVector(SequencesDatabase, OutPut, Label);
		//cout << "total pattern vector size is " << TotalPatternVector.size() << endl;
		JensenShannonDivergence(TotalPatternVector, Label, outputDM + "/DM.txt");
		//CosineSimilarityMatirx(TotalPatternVector, Label, DMPath + "/6CLIDFPVCSM.txt");


		//TotalPatternVector = IDFWeightPatternVector(SequencesDatabase, TempOutPut, Label, DMPath + "/TFIDFPV.txt");
		//JensenShannonDivergence(TotalPatternVector, Label, DMPath + "/6TFIDFPVJSDM.txt");
		//CosineSimilarityMatirx(TotalPatternVector, Label, DMPath + "/8TFIDFPVCSM.txt");
	}

	cout << "output  " << OutPut.size() << endl;
	cout << "tempoutput  " << TempOutPut.size() << endl;
	return 0;
}


//using namespace std;
// parameters are projected database, original input database,
// support ratio, minimal non-wildcard items,
// maximal length of continuous wildcard, where to output results
void mining(ProjectedDatabase & TempProData,
	const vector <string>& SequencesDatabase,
	const double & MinSupRatio,
	const int & Min_Pat_Length,
	const int & Max_WildCard_Length,
	int & MaxPL,
	vector <MotifSupportPosition> & TempOutPut,
	vector <MotifSupportPosition> & OutPut)
{
	if (TempProData.GetSupport() < SequencesDatabase.size() * MinSupRatio)
		return;  // end current projected database if support is low
				 // if support is more than MinSupRatio and non-wildcard items is
				 // more than Min_Pat_Length, output pattern
	//cout << "mining 1" << endl;

	if ((TempProData.GetPrefixSize() - TempProData.GetTotalWildCardLength()) >= Min_Pat_Length)
	{

		/*for (int i = 0; i <= TempProData.GetPrefixSize() - 1; i++)
		{
			cout << TempProData.Prefix[i];
		}
		cout << endl;*/

		MotifSupportPosition TempMotifSupportPosition;
		TempMotifSupportPosition.Motif = TempProData.GetPrefix();
//		cout << "mining 2" << endl;
		if (TempMotifSupportPosition.Motif[TempMotifSupportPosition.Motif.size() - 1] != 'x') 
		{

			/*for (int i = 0; i <= TempProData.GetPrefixSize() - 1; i++)
			{
				cout << TempProData.Prefix[i];
			}
			cout << endl;*/
			TempMotifSupportPosition.Support = TempProData.GetSupport();
			//kyy************************************************************************************************
			for (int proteinindex = 0; proteinindex <= TempProData.ForEachProtein.size() - 1; proteinindex++)
			{
				RPIPP TempRPIPP;
				TempRPIPP.RPI = TempProData.ForEachProtein[proteinindex].RealProteinIndex;
				for (int positionindex = 0; positionindex <= TempProData.ForEachProtein[proteinindex].CurrentPrefix.size() - 1; positionindex++)
				{
					TempRPIPP.PP.push_back(TempProData.ForEachProtein[proteinindex].CurrentPrefix[positionindex]);
				}
				TempMotifSupportPosition.OUTRPIPP.push_back(TempRPIPP);
			}
			//kyy************************************************************************************************
			TempOutPut.push_back(TempMotifSupportPosition);

			if (TempOutPut.size() > 2)
			{
				if (TempOutPut[TempOutPut.size() - 1].Support == TempOutPut[TempOutPut.size() - 2].Support)
				{
					if (TempOutPut[TempOutPut.size() - 1].Motif[0] != TempOutPut[TempOutPut.size() - 2].Motif[0]
						|| TempOutPut[TempOutPut.size() - 1].Motif.size() <= TempOutPut[TempOutPut.size() - 2].Motif.size())
					{
						OutPut.push_back(TempOutPut[TempOutPut.size() - 2]);
					}
				}
				else
				{
					OutPut.push_back(TempOutPut[TempOutPut.size() - 2]);
				}
				if (OutPut[OutPut.size() - 1].Motif.size() > MaxPL)
				{
					MaxPL = OutPut[OutPut.size() - 1].Motif.size();
				}
				//cout << "maxpl is " << MaxPL << endl;
			}
		}
	}
//	cout << "mining 3" << endl;
	if (TempProData.GetCurrentWildCardLength() < Max_WildCard_Length)
	{
		/*for (int i = 0; i <= TempProData.GetPrefixSize() - 1; i++)
		{
			cout << TempProData.Prefix[i];
		}
		cout << endl;*/


		ProjectedDatabase ProData(TempProData, false);
		ProData.MoveAll();
		mining(ProData, SequencesDatabase, MinSupRatio, Min_Pat_Length, Max_WildCard_Length, MaxPL, TempOutPut, OutPut);
	}
//	cout << "mining 4" << endl;
	for (int AAIndex = 19; AAIndex >= 0; AAIndex--)
	{
		/*for (int i = 0; i <= TempProData.GetPrefixSize() - 1; i++)
		{
			cout << TempProData.Prefix[i];
		}
		cout << endl;*/


		ProjectedDatabase ProData(TempProData, true);
		ProData.UpdateProData(TempProData, AASubType[AAIndex], SequencesDatabase);
		mining(ProData, SequencesDatabase, MinSupRatio, Min_Pat_Length, Max_WildCard_Length, MaxPL, TempOutPut, OutPut);
	}
}


vector<vector<int> > BinaryPatternVector( const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label, string OutFileName)
{

	int PatternNumbers = OutPut.size();
	vector<vector<int> > RefPatternVector;
	vector<int> TempRefPatternVector(PatternNumbers, 0);
	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
		RefPatternVector.push_back(TempRefPatternVector);
	}

	for (int PatternIndex = 0; PatternIndex <= OutPut.size() - 1; PatternIndex++)
	{//for each pattern
		for (int RPIIndex = 0; RPIIndex <= OutPut[PatternIndex].OUTRPIPP.size() - 1; RPIIndex++)
		{ //for each RPI
			RefPatternVector[OutPut[PatternIndex].OUTRPIPP[RPIIndex].RPI][PatternIndex]
				= 1;
		}
	}

	return RefPatternVector;
}


vector<vector<int> > IDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label, string OutFileName)
{

	int PatternNumbers = OutPut.size();
	vector<vector<int> > RefPatternVector;
	vector<int> TempRefPatternVector(PatternNumbers, 0);
	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
		RefPatternVector.push_back(TempRefPatternVector);
	}

	for (int PatternIndex = 0; PatternIndex <= OutPut.size() - 1; PatternIndex++)
	{//for each pattern
		for (int RPIIndex = 0; RPIIndex <= OutPut[PatternIndex].OUTRPIPP.size() - 1; RPIIndex++)
		{ //for each RPI
			RefPatternVector[OutPut[PatternIndex].OUTRPIPP[RPIIndex].RPI][PatternIndex]
				= log2(SequencesDatabase.size() / OutPut[PatternIndex].Support);
		}
	}

	return RefPatternVector;
}

vector<vector<int> > FIDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label, string OutFileName)
{
	int PatternNumbers = OutPut.size();
	vector<vector<int> > RefPatternVector;
	vector<int> TempRefPatternVector(PatternNumbers, 0);
	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
		RefPatternVector.push_back(TempRefPatternVector);
	}

	for (int PatternIndex = 0; PatternIndex <= OutPut.size() - 1; PatternIndex++)
	{//for each pattern
		for (int RPIIndex = 0; RPIIndex <= OutPut[PatternIndex].OUTRPIPP.size() - 1; RPIIndex++)
		{ //for each RPI
			RefPatternVector[OutPut[PatternIndex].OUTRPIPP[RPIIndex].RPI][PatternIndex]
				= (MaxPatternLength / OutPut[PatternIndex].Motif.size()) * log2(SequencesDatabase.size() / OutPut[PatternIndex].Support);
		}
	}

	return RefPatternVector;
}



vector<vector<int> > WIDFWeightPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label)
{

	int PatternNumbers = OutPut.size();
	vector<vector<int> > RefPatternVector;
	vector<int> TempRefPatternVector(PatternNumbers, 0);
	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
		RefPatternVector.push_back(TempRefPatternVector);
	}

	for (int PatternIndex = 0; PatternIndex <= OutPut.size() - 1; PatternIndex++)
	{//for each pattern
		for (int RPIIndex = 0; RPIIndex <= OutPut[PatternIndex].OUTRPIPP.size() - 1; RPIIndex++)
		{ //for each RPI
			RefPatternVector[OutPut[PatternIndex].OUTRPIPP[RPIIndex].RPI][PatternIndex]
				= (OutPut[PatternIndex].Motif.size()) * log2(SequencesDatabase.size() / OutPut[PatternIndex].Support);
		}
	}

	return RefPatternVector;
}



vector<vector<int> > ConsiderWeightDuplicateOverlapPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label, string OutFileName)
{


	vector<vector<int> > ReferenceSequenceDatabase;
	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
		vector<int> TempPerSequenceDatabase(SequencesDatabase[SequenceDatabaseIndex].size(), 0);// 2017-11-03 change
																								//vector<int> TempPerSequenceDatabase((SequencesDatabase[SequenceDatabaseIndex].size() - 1), 0);
		ReferenceSequenceDatabase.push_back(TempPerSequenceDatabase);
	}
	//kyy*****************************************************************
	/*for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
	{
	cout << SequencesDatabase[SequenceDatabaseIndex].size() - 1 << endl;
	cout << ReferenceSequenceDatabase[SequenceDatabaseIndex].size() << endl;
	for (int i = 0; i <= ReferenceSequenceDatabase[SequenceDatabaseIndex].size() - 1; i++)
	{
	cout << ReferenceSequenceDatabase[SequenceDatabaseIndex][i];
	}
	cout << endl;
	}*/
	//kyy*****************************************************************

	vector<int> OnePatternVector;
	vector<vector<int> > TotalPatternVector;

	for (int SequenceIndex = 0; SequenceIndex <= SequencesDatabase.size() - 1; SequenceIndex++)
	{//for each sequences
		for (int OutPutPatternIndex = 0; OutPutPatternIndex <= OutPut.size() - 1; OutPutPatternIndex++)
		{ //for each pattern 
			int PatternOccur = 0;
			int OccurCount = 0;
			for (int OUTRPIPPIndex = 0; OUTRPIPPIndex <= OutPut[OutPutPatternIndex].OUTRPIPP.size() - 1; OUTRPIPPIndex++)// for each RPI-PP ( times one pattern occur in one suquence for all suquences with this pattern ) 
			{
				if (SequenceIndex == OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].RPI)//pattern occur in this suquence 
				{
					PatternOccur++;

					for (int PPIndex = 0; PPIndex <= OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].PP.size() - 1; PPIndex++)//for each pattern occur position 
					{
						int Number1Count = 0;
						for (int MarkIndex = (OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].PP[PPIndex] - OutPut[OutPutPatternIndex].Motif.size() + 1);
							MarkIndex <= OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].PP[PPIndex];
							MarkIndex++)
						{
							/*for (int MarkIndex = OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].PP[PPIndex];
							MarkIndex >= (OutPut[OutPutPatternIndex].OUTRPIPP[OUTRPIPPIndex].PP[PPIndex] - OutPut[OutPutPatternIndex].Motif.size() + 1);
							MarkIndex--)
							{*/
							if (ReferenceSequenceDatabase[SequenceIndex][MarkIndex] == 0)
							{
								ReferenceSequenceDatabase[SequenceIndex][MarkIndex] = 1;
								Number1Count++;
							}
						}
						OccurCount = OccurCount + Number1Count;
					}
				}
			}
			if (PatternOccur != 0)
				OnePatternVector.push_back(OccurCount);
			else OnePatternVector.push_back(0);
		}
		TotalPatternVector.push_back(OnePatternVector);
		OnePatternVector.clear();
	}
	//	ofstream ConsiderWeightDuplicateOverlapPatternVector(OutFileName.c_str());
	//	for (int i = 0; i <= TotalPatternVector.size() - 1; i++)
	//	{
	//		ConsiderWeightDuplicateOverlapPatternVector << Label[i] << endl;
	//		for (int j = 0; j <= TotalPatternVector[i].size() - 1; j++)
	//		{
	//			ConsiderWeightDuplicateOverlapPatternVector << TotalPatternVector[i][j] << " ";
	//		}
	//		ConsiderWeightDuplicateOverlapPatternVector << endl;
	//	}
	return TotalPatternVector;
}


//vector<vector<int> > ConsiderWeightDuplicateOverlapPatternVector(const vector<string> &SequencesDatabase, const vector<MotifSupportPosition> &OutPut, const vector<string> &Label, string OutFileName)
//{
//
//	vector<vector<int> > ReferenceSequenceDatabase;
//	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
//	{
//		vector<int> TempPerSequenceDatabase(SequencesDatabase[SequenceDatabaseIndex].size(), 0);																								
//		ReferenceSequenceDatabase.push_back(TempPerSequenceDatabase);
//	}
//
//	int PatternNumbers = OutPut.size();
//	vector<vector<int> > RefPatternVector;
//	vector<int> TempRefPatternVector(PatternNumbers, 0);
//	for (int SequenceDatabaseIndex = 0; SequenceDatabaseIndex <= SequencesDatabase.size() - 1; SequenceDatabaseIndex++)
//	{
//		RefPatternVector.push_back(TempRefPatternVector);
//	}
//	 
//	for (int PatternIndex = 0; PatternIndex <= OutPut.size() - 1; PatternIndex++)
//	{//for each pattern
//		for (int RPIIndex = 0; RPIIndex <= OutPut[PatternIndex].OUTRPIPP.size() - 1; RPIIndex++)
//		{ //for each RPI
//			
//			int OccurCount = 0;
//			for (int PPIndex = 0; PPIndex <= OutPut[PatternIndex].OUTRPIPP[RPIIndex].PP.size() - 1; PPIndex++)//for each pattern occur position 
//			{
//				int Number1Count = 0;
//				for (int MarkIndex = (OutPut[PatternIndex].OUTRPIPP[RPIIndex].PP[PPIndex] - OutPut[PatternIndex].Motif.size() + 1);
//					MarkIndex <= OutPut[PatternIndex].OUTRPIPP[RPIIndex].PP[PPIndex]; MarkIndex++)
//				{
//					if (ReferenceSequenceDatabase[RPIIndex][MarkIndex] == 0)
//					{
//						ReferenceSequenceDatabase[RPIIndex][MarkIndex] = 1;
//						Number1Count++;
//					}
//				}
//				OccurCount = OccurCount + Number1Count;
//				cout << "occurcount is " << OccurCount << endl;
//			}	
//			
//			RefPatternVector[OutPut[PatternIndex].OUTRPIPP[RPIIndex].RPI][PatternIndex] = OccurCount;
//		}
//	}
//
//	return RefPatternVector;
//}




void JaccardDistanceMatrix( const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName)
{
	vector< vector<double> > Matrix;
	vector< double > PairDistance;
	double JaccardDistance;
	for (int FirstBinarySequenceIndex = 0; FirstBinarySequenceIndex < TotalPatternVector.size(); FirstBinarySequenceIndex++)//the first sequence as reference sequence;
	{
		for (int NextBinarySequenceIndex = 0; NextBinarySequenceIndex < TotalPatternVector.size(); NextBinarySequenceIndex++)//another sequence compared with the first sequence
		{
			double a = 0; double b = 0; double c = 0;
			for (int BinarySequenceBitIndex = 0; BinarySequenceBitIndex <= TotalPatternVector[0].size() - 1; BinarySequenceBitIndex++)//compare each bit
			{
				if (TotalPatternVector[FirstBinarySequenceIndex][BinarySequenceBitIndex] == 1 && TotalPatternVector[NextBinarySequenceIndex][BinarySequenceBitIndex] == 0)
					a++;
				if (TotalPatternVector[FirstBinarySequenceIndex][BinarySequenceBitIndex] == 0 && TotalPatternVector[NextBinarySequenceIndex][BinarySequenceBitIndex] == 1)
					b++;
				if (TotalPatternVector[FirstBinarySequenceIndex][BinarySequenceBitIndex] == 1 && TotalPatternVector[NextBinarySequenceIndex][BinarySequenceBitIndex] == 1)
					c++;
			}
			if (a + b + c == 0)
				JaccardDistance = 0;
			else JaccardDistance = (a + b) / (a + b + c);//1008 this is the distance of two sequence
			//else JaccardDistance = c / (a + b + c);//1008 this is the similarity of two sequence
			PairDistance.push_back(JaccardDistance);
			a = 0; b = 0; c = 0;
		}
		Matrix.push_back(PairDistance);
		PairDistance.clear();
	}
	ofstream BinaryPatternVectorJaccardDistanceMatrix(OutFileName.c_str());
	BinaryPatternVectorJaccardDistanceMatrix << Label.size() << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		BinaryPatternVectorJaccardDistanceMatrix << Label[i] << "     ";
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			BinaryPatternVectorJaccardDistanceMatrix << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << Matrix[i][j];
			BinaryPatternVectorJaccardDistanceMatrix << " ";
		}
		BinaryPatternVectorJaccardDistanceMatrix << endl;
	}
}


void TanimotoDistanceMatrix(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, string OutFileName)
{
	vector<vector<double> >  Matrix;
	vector<double> PairDistance;
	double distance;
	int PatternVectorLength = TotalPatternVector[0].size();
	for (int RefPatternVectorIndex = 0; RefPatternVectorIndex <= TotalPatternVector.size() - 1; RefPatternVectorIndex++)
	{
		for (int AnotherPatternVectorIndex = 0; AnotherPatternVectorIndex <= TotalPatternVector.size() - 1; AnotherPatternVectorIndex++)
		{
			double Numerator = 0;
			double Denominator = 0;
			double RefSquare = 0;
			double AnotherSquare = 0;
			for (int BitIndex = 0; BitIndex <= PatternVectorLength - 1; BitIndex++)
			{
				Numerator = Numerator + TotalPatternVector[RefPatternVectorIndex][BitIndex] * TotalPatternVector[AnotherPatternVectorIndex][BitIndex];
				RefSquare = RefSquare + pow(TotalPatternVector[RefPatternVectorIndex][BitIndex], 2);
				AnotherSquare = AnotherSquare + pow(TotalPatternVector[AnotherPatternVectorIndex][BitIndex], 2);
			}
			Denominator = RefSquare + AnotherSquare - Numerator;
			if (Denominator == 0 )
				distance = 0;
			else distance = 1 - Numerator / Denominator;//1008 is it a distance?
			//distance = Numerator / Denominator; // 1008 Tanimoto coefficient or Tanimoto distance
			PairDistance.push_back(distance);
		}
		Matrix.push_back(PairDistance);
		PairDistance.clear();
	}
	ofstream BinaryPatternVectorTanimotoDistanceMatrix(OutFileName.c_str());
	BinaryPatternVectorTanimotoDistanceMatrix << Label.size() << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		BinaryPatternVectorTanimotoDistanceMatrix << Label[i] << "     ";
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			BinaryPatternVectorTanimotoDistanceMatrix << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << Matrix[i][j];
			BinaryPatternVectorTanimotoDistanceMatrix << " ";
		}
		BinaryPatternVectorTanimotoDistanceMatrix << endl;
	}
}

void JensenShannonDivergence(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label, const string OutFileName)
{
	vector< vector<double> > Matrix;
	vector< double > PairDistance;
	double P1, P2, JSD_1 = 0, JSD_2 = 0, average = 0, JSD = 0;
	
	for (int RefPatternVectorIndex = 0; RefPatternVectorIndex <= TotalPatternVector.size() - 1; RefPatternVectorIndex++)
	{
		double RefSum = accumulate(TotalPatternVector[RefPatternVectorIndex].begin(), TotalPatternVector[RefPatternVectorIndex].end(), 0);
		//cout << "refsum " << RefSum << endl;
		/*if (RefSum == 0){
			JSD = 1;
			Matrix.push_back(TotalPatternVector[RefPatternVectorIndex].size(), 1);
			Matrix[RefPatternVectorIndex][0] = 0;
			continue;
		}*/
		for (int AnotherPatternVectorIndex = 0; AnotherPatternVectorIndex <= TotalPatternVector.size() - 1; AnotherPatternVectorIndex++)
		{
			double AnoSum = accumulate(TotalPatternVector[AnotherPatternVectorIndex].begin(), TotalPatternVector[AnotherPatternVectorIndex].end(), 0);
			//cout << "anosum " << AnoSum << endl;
			if (RefPatternVectorIndex == AnotherPatternVectorIndex)
			{
				PairDistance.push_back(JSD);
			}
			else if (RefSum == 0 || AnoSum == 0)
			{
				JSD = 1;
				PairDistance.push_back(JSD);
				JSD = 0;
			}
			else {
				for (int BitIndex = 0; BitIndex <= TotalPatternVector[0].size() - 1; BitIndex++)
				{
					P1 = TotalPatternVector[RefPatternVectorIndex][BitIndex] / RefSum;
					//cout << "p1 " << P1 << endl;
					P2 = TotalPatternVector[AnotherPatternVectorIndex][BitIndex] / AnoSum;
					//cout << "p2 " << P2 << endl;
					average = (P1 + P2) / 2;
					if (!P1 && !P2) {
						continue;
					}
					else if (!P1) {
						JSD_2 -= P2;
						continue;
					}
					else if (!P2) {
						JSD_2 -= P1;
						continue;
					}
					else {
						JSD_1 += -P1 * log2(average / P1);
						JSD_2 += -P2 * log2(average / P2);
					}
					//cout << "jsd1 " << JSD_1 << endl;
					//cout << "jsd2 " << JSD_2 << endl;
				}
				JSD = fabs(0.5*JSD_1 + 0.5*JSD_2);
				//cout << "jsd " << JSD << endl;
				PairDistance.push_back(JSD);
				JSD_1 = 0;
				JSD_2 = 0;
				JSD = 0;
			}					
		}	
			Matrix.push_back(PairDistance);
			PairDistance.clear();	
	}

	ofstream PatternVectorJSDMatrix(OutFileName.c_str());
	PatternVectorJSDMatrix << Label.size() << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		PatternVectorJSDMatrix << Label[i] << "     ";
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			PatternVectorJSDMatrix << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << Matrix[i][j];
			PatternVectorJSDMatrix << " ";
		}
		PatternVectorJSDMatrix << endl;
	}

}





void EuclideanDistanceMatrix(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label,const string OutFileName)
{
	vector<vector<double> >  Matrix;
	vector<double> PairDistance;
	double distance;
	int PatternVectorLength = TotalPatternVector[0].size();
	for (int RefPatternVectorIndex = 0; RefPatternVectorIndex <= TotalPatternVector.size() - 1; RefPatternVectorIndex++)
	{
		for (int AnotherPatternVectorIndex = 0; AnotherPatternVectorIndex <= TotalPatternVector.size() - 1; AnotherPatternVectorIndex++)
		{
			int sum = 0;
			for (int BitIndex = 0; BitIndex <= PatternVectorLength - 1; BitIndex++)
			{
				sum = sum + pow((TotalPatternVector[RefPatternVectorIndex][BitIndex] - TotalPatternVector[AnotherPatternVectorIndex][BitIndex]), 2);
			}
			distance = sqrt(sum);//1008 this is the distance
			//if(sqrt(sum)==0) distance = 0;// 1008 this is the distance-similarity
			//else distance = 1 / sqrt(sum);// 1008 this is the distance-similarity
			PairDistance.push_back(distance);
		}
		Matrix.push_back(PairDistance);
		PairDistance.clear();
	}
	ofstream PatternVectorEuclideanDistanceMatrix(OutFileName.c_str());
	PatternVectorEuclideanDistanceMatrix << "   " << Label.size() << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		PatternVectorEuclideanDistanceMatrix << Label[i] << "     ";
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			PatternVectorEuclideanDistanceMatrix << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << Matrix[i][j];
			PatternVectorEuclideanDistanceMatrix << " ";
		}
		PatternVectorEuclideanDistanceMatrix << endl;
	}
}

void CosineSimilarityMatirx(const vector<vector<int> > &TotalPatternVector, const vector<string> &Label,const string OutFileName)
{
	vector<vector<double> >  Matrix;
	vector<double> PairDistance;
	double distance;
	int PatternVectorLength = TotalPatternVector[0].size();
	for (int RefPatternVectorIndex = 0; RefPatternVectorIndex <= TotalPatternVector.size() - 1; RefPatternVectorIndex++)
	{
		for (int AnotherPatternVectorIndex = 0; AnotherPatternVectorIndex <= TotalPatternVector.size() - 1; AnotherPatternVectorIndex++)
		{
			double Numerator = 0;
			double Denominator = 0;
			double RefSquare = 0;
			double AnotherSquare = 0;
			double RefEuclideanNorm = 0;
			double AnotherEuclideanNorm = 0;
			for (int BitIndex = 0; BitIndex <= PatternVectorLength - 1; BitIndex++)
			{
				Numerator = Numerator + (TotalPatternVector[RefPatternVectorIndex][BitIndex] * TotalPatternVector[AnotherPatternVectorIndex][BitIndex]);
				//cout << Numerator << endl;
				RefSquare = RefSquare + pow(TotalPatternVector[RefPatternVectorIndex][BitIndex], 2);
				AnotherSquare = AnotherSquare + pow(TotalPatternVector[AnotherPatternVectorIndex][BitIndex], 2);
			}
			RefEuclideanNorm = sqrt(RefSquare);
			AnotherEuclideanNorm = sqrt(AnotherSquare);
			Denominator = RefEuclideanNorm * AnotherEuclideanNorm;
			if (Denominator==0)
				distance=0;
			else distance = fabs(1 - Numerator / Denominator);//1008 is it a distance?
			//distance = Numerator / Denominator;// 1008 this is the cosine similarity
			PairDistance.push_back(distance);
		}
		Matrix.push_back(PairDistance);
		PairDistance.clear();
	}
	ofstream PatternVectorCosineSimilarityDistanceMatrix(OutFileName.c_str());
	PatternVectorCosineSimilarityDistanceMatrix << Label.size() << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		PatternVectorCosineSimilarityDistanceMatrix << Label[i] << "     ";
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			PatternVectorCosineSimilarityDistanceMatrix << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << Matrix[i][j];
			PatternVectorCosineSimilarityDistanceMatrix << " ";
		}
		PatternVectorCosineSimilarityDistanceMatrix << endl;
	}
}










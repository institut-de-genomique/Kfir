#include <iostream>
#include <fstream>
#include <cstdio>

#include <vector>
#include <algorithm>
#include <utility>
#include <time.h>
#include <malloc.h>

#include "dust.h"
#include "DnaDictionary.h"
#include "ReadFile.h"
#include <gzstream.h>
#include <libgen.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <math.h>


using namespace std;
using namespace boost::algorithm;

void usage();
void error (string msg);
bool testLowComplexity(string read, s32 testComplexity, f8 threshold);
bool mask_read(string read);
f8 loga2(float nb);
bool shanEnt(string read, f8 threshold);
bool check_read(DnaDictionary &dict, string read, s32 nb_words);
bool existence(char* name);


int main(int argc, char** argv) 
{

  s32 c, word_size=19, nb_words=1;
  f8 threshold = 0;
  char *genomeFilename = NULL, *readFilename1 = NULL, *readFilename2 = NULL, *outputPrefix = NULL ;
  s32 testComplexity = 0; //takes values : 1 if dust, 2 if shannon entropy
  bool testComplexityFlag = 0;
  string outputString;
  
  while((c = getopt(argc, argv, "g:r:p:k:n:o:dt::h")) != -1) {
    switch (c) {
    case 'g':
      genomeFilename = optarg;
      break;
    case 'r':
      readFilename1 = optarg;
      break;
    case 'p':
      readFilename2 = optarg;
      break;
    case 'k':
      word_size = atoi(optarg);
      break;
    case 'n':
      nb_words = atoi(optarg);
      break;
    case 'o':
      outputPrefix = optarg;
      break;
    case 'd':
      if (testComplexityFlag) error("Need to specify only one complexity test : Dust OR Shannon Entropy.");
      testComplexityFlag = 1;
      testComplexity = 1;
      break;
    case 't':
      if (testComplexityFlag) error("Need to specify only one complexity test : Dust OR Shannon Entropy.");
      testComplexityFlag = 1;
      if (optarg == NULL) { threshold = 1.59; testComplexity = 2; }
      else { threshold = atof(optarg); testComplexity = 2; }
      break;
    case 'h': 
      usage();
      break;
    default :
      abort();
  }
  }
 
  if (genomeFilename == NULL) error("Need to provide a reference genome for the dictionary.");
  else if (!existence(genomeFilename)) error("Need to provide a valid filename for reference genome.");
  if (readFilename1 == NULL) error("Need to provide a reads filename.");
  if (readFilename2 != NULL && readFilename1 == readFilename2) error("Need to provide two distinct reads filenames.");
  if (outputPrefix == NULL) outputString = "";
  else { outputString = (string)outputPrefix; }


  string strReadFilename1 = string(readFilename1), strReadFilename2;
  vector<string> listOfReads1, listOfReads2;
  split(listOfReads1, strReadFilename1, is_any_of(",;"));
  if (readFilename2 != NULL) {
    strReadFilename2 = string(readFilename2);
    split(listOfReads2, strReadFilename2, is_any_of(",;"));
    
    if (listOfReads1.size() != listOfReads2.size()) error("Need to provide the same number of files for read1 and read2");
  }  

  // CREATION OF A WORD DICTIONARY - NEED TO SPECIFY THE SIZE OF WORDS

  DnaDictionary dict(word_size);
  ReadFile genome(genomeFilename); 
  genome.loadAndCount(dict); 
  cerr<<"Genome loaded successfully. Total number of different words :  "<< dict.getNbDiffWords() << endl;
 


  // FILTERING ON A SINGLE READ FILE

  if (!readFilename2) {

    vector<string>::iterator itList;
    s32 processedReads = 0;
    for (itList = listOfReads1.begin(); itList != listOfReads1.end(); itList++){
      
      //open the file, create a corresponding output file, and read through it
      igzstream _fin1; 
      fstream _fout1;
      string read1_S1, read1_S2, read1_S3, read1_S4; // read identifier, read sequence, interlign, read quality
      
      char* charItList = strdup((*itList).c_str());
      if (!existence(charItList)) {
	cerr << "[Warning] File not found : "+(string)charItList << endl;
	continue;}
      _fin1.open(charItList, ios_base::in | ios_base::binary);
  
      string read1_filtered = outputString +  string( basename(charItList) );
      u8 pos = read1_filtered.find(".gz");
      read1_filtered = read1_filtered.substr(0,pos);
      
      free(charItList);
      _fout1.open(read1_filtered.c_str(), ios_base::out | ios_base::binary);

      _fin1>>read1_S1; 
      _fin1>>read1_S2; 
      _fin1>>read1_S3;
      _fin1>>read1_S4; 
      
      while (!_fin1.eof()) {
	
	if (testComplexity && check_read(dict, read1_S2, nb_words)){
	  if (testLowComplexity(read1_S2, testComplexity, threshold)) {
	    _fout1<<read1_S1<<endl;
	    _fout1<<read1_S2<<endl;
	    _fout1<<read1_S3<<endl;
	    _fout1<<read1_S4<<endl;
	  }
	}
	else if (check_read(dict, read1_S2, nb_words)){
	  _fout1<<read1_S1<<endl;
	  _fout1<<read1_S2<<endl;
	  _fout1<<read1_S3<<endl;
	  _fout1<<read1_S4<<endl;
	}
	_fin1>>read1_S1; 
	_fin1>>read1_S2;  
	_fin1>>read1_S3;
	_fin1>>read1_S4; 
      }
      _fin1.close();
      _fout1.close();
      processedReads++;
    }
    cerr << "A total of " << processedReads << " reads file(s) were processed."<< endl;
    return 0;
  }


  // FILTERING ON PAIRED-END READ FILES

  else {    
    s32 processedReads = 0;
    for(s8 itList=0; itList<(s8)listOfReads1.size(); itList++){
      
      igzstream _fin1, _fin2;
      fstream _fout1, _fout2;
      string read1_S1, read1_S2, read1_S3, read1_S4, read2_S1, read2_S2, read2_S3, read2_S4;
      char *charItList1 = strdup(listOfReads1[itList].c_str()), *charItList2 = strdup(listOfReads2[itList].c_str());
      if (!existence(charItList1)) {
	cerr << "[Warning] File not found : "+(string)charItList1 << endl;
	continue;}
      if (!existence(charItList2)) {
	cerr << "[Warning] File not found : "+(string)charItList2 << endl;
	continue;}
      _fin1.open(charItList1, ios_base::in | ios_base::binary);
      _fin2.open(charItList2, ios_base::in | ios_base::binary);
      
      string read1_filtered = outputString + (string)basename(charItList1), read2_filtered = outputString+(string)basename(charItList2);
     
      int pos1 = read1_filtered.find(".gz"), pos2 = read2_filtered.find(".gz");
      read1_filtered = read1_filtered.substr(0, pos1);
      read2_filtered = read2_filtered.substr(0, pos2);
      

      _fout1.open(read1_filtered.c_str(), ios_base::out | ios_base::binary);
      _fout2.open(read2_filtered.c_str(), ios_base::out | ios_base::binary);
      
      
      free(charItList1);
      free(charItList2);

      _fin1>>read1_S1;
      _fin1>>read1_S2;
      _fin1>>read1_S3;
      _fin1>>read1_S4;
      
      _fin2>>read2_S1;
      _fin2>>read2_S2;
      _fin2>>read2_S3;
      _fin2>>read2_S4;

      while ((!_fin1.eof()) || (!_fin2.eof())) { //redondant car les 2 fichiers sont synchronises et doivent faire la meme taille
	
	if (testComplexity && (check_read(dict, read1_S2, nb_words) || check_read(dict, read2_S2, nb_words) )) {

	  if (testLowComplexity(read1_S2, testComplexity, threshold) || testLowComplexity(read2_S2, testComplexity, threshold)) {
	    _fout1<<read1_S1<<endl;
	    _fout1<<read1_S2<<endl;
	    _fout1<<read1_S3<<endl;
	    _fout1<<read1_S4<<endl;
	    
	    _fout2<<read2_S1<<endl;
	    _fout2<<read2_S2<<endl;
	    _fout2<<read2_S3<<endl;
	    _fout2<<read2_S4<<endl;
	  }
	}

	else if (check_read(dict, read1_S2, nb_words) || check_read(dict, read2_S2, nb_words)) {
	  _fout1<<read1_S1<<endl;
	  _fout1<<read1_S2<<endl;
	  _fout1<<read1_S3<<endl;
	  _fout1<<read1_S4<<endl;
	  
	  _fout2<<read2_S1<<endl;
	  _fout2<<read2_S2<<endl;
	  _fout2<<read2_S3<<endl;
	  _fout2<<read2_S4<<endl;
	}	
	_fin1>>read1_S1;
	_fin1>>read1_S2;
	_fin1>>read1_S3;
	_fin1>>read1_S4;
	
	_fin2>>read2_S1;
	_fin2>>read2_S2;
	_fin2>>read2_S3;
	_fin2>>read2_S4;

      }
      _fin1.close();
      _fin2.close();
      _fout1.close();
      _fout2.close();
      
      processedReads++;
    }
    cerr << "A total of "<< processedReads << " reads files pairs were processed." << endl;
    return 0;
  }
}


// to test the low complexity (dust or shannon entropy)
bool testLowComplexity(string read, s32 testComplexity, f8 threshold = NULL) {
  if (testComplexity == 1) return mask_read(read);
  else if (testComplexity == 2) return shanEnt(read, threshold);
  else return 1;
}

// to mask the low complexity segments of a read
bool mask_read(string read) {
  s32 ratio = 75;
  s32 unmasked_len = 30;;
  char* read_chr = strdup(read.c_str());
  s32 size = read.size();
  s32 sum = 0;
  s32 *psum = &sum;
  dust(size, read_chr, psum);
  float frac = ((float)sum/size)*100;
  free(read_chr);
  if ((frac < ratio) && ((size-sum) > unmasked_len)) {return TRUE;}
  else {return FALSE;}
}


// to calculate the logarithm
f8 loga2(float nb){
  if (nb==0) {return 0;}
  else {return log2(nb);}
}

// to calculate the Shannon entropy
bool shanEnt(string read, f8 threshold){
  f8 entropy;
  f8 propA=0, propC=0, propG=0, propT=0;
  for (s32 i = 0; i < (s32)read.size(); i++){
    switch( tolower(read[i]) )
      {
      case 'a' : 
	propA++;
	break;
      case 'c' : 
	propC++;
	break;
      case 'g' : 
	propG++;
	break;
      case 't' : 
	propT++;
	break;
      default :
	break;
      }
  }
  propA/=read.size();
  propC/=read.size();
  propG/=read.size();
  propT/=read.size();
  
  entropy = -(propA * loga2(propA) + propC * loga2(propC) + propG * loga2(propG) + propT * loga2(propT));
  
  if (entropy < threshold) {return FALSE;}
  else{ return TRUE;}
}

// to verify if a read belongs to a genome
bool check_read(DnaDictionary &dict, string read, s32 nb_words){
  s32 count_words = 0;
  if ((s32)read.size() < dict.getWordSize()) return FALSE; 
  for (s32 i = 0; i < (s32)read.size() - dict.getWordSize() + 1; i++) {
    string word = read.substr(i, dict.getWordSize()); 
    if((dict.existWord(word)) && (++count_words >= nb_words)) {
      return TRUE;}
  }
  return FALSE;
}


// to send an error message
void error(string msg) {
  cerr << "[Error] " << msg << endl;
  cerr << "See filter -h for more details." << endl;
  exit(1);
}

// to test the existence of a file
bool existence(char* name) {
  ifstream f(name);
  if (f.good()) {
    f.close();
    return true;
  }
  else {
    f.close();
    return false;
  }   
}

// the result of kfir -h
void usage() {
	cerr << "--------------------------------------------------------------------------------------------" << endl;
	cerr << "Compares a set of reads to a genome and yields the subset belonging to it" << endl << endl;
	cerr << "--------------------------------------------------------------------------------------------" << endl;
	cerr << "Usage : kfir -g <genome file> -r <reads file> -o <output prefix> {Options}" << endl;
	cerr << "!! Note : Arguments with * are required, the other are optionnal !!" << endl << endl;
	cerr << "  INPUT FILES" << endl;
	cerr << "*     -g  <str> : reference genome/chromosome/scaffold" << endl;
	cerr << "*     -r  <str> : first file of reads - or list of files (.fq or .fq.gz)" << endl;
	cerr << "      -p  <str> : second file of reads if paired-end files - or list of files (.fq or .fq.gz)" << endl;
	cerr << "      -k  <int> : word size for the dictionary, default is 19" << endl;
	cerr << "      -n  <int> : number of words in common to validate a read as belonging to the genome, default is 1" << endl;
	cerr << "      -o  <str> : prefix of the output files" << endl;
	cerr << "  LOW COMPLEXITY TEST : choose one of the following options, default no complexity test" << endl;
	cerr << "      -d        : add a dust complexity mask on the read files" << endl;
	cerr << "      -t<float> : if missing, no Shannon entropy comparison, if -t default threshold is 1.59" << endl;
	cerr << "      -h        : this help" << endl;
	cerr << "--------------------------------------------------------------------------------------------" << endl;
	exit(1);
}

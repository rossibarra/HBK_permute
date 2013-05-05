//Hudson, Boos, and Kaplan (1992) Mol. Biol. Evol. 9(1):138-151
//Fst is calculated according to equation 3a of
//Charlesworth (1998) Mol. Biol. Evol.: 15:538-543
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <algorithm>
#include <numeric>
#include <getopt.h>

#if defined(__GNUG__ ) &&  __GNUC__ >= 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/PolySites.hpp>
#include <Sequence/FST.hpp>
#include <RandomNumbers.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/stateCounter.hpp>

using namespace Sequence;
using namespace Alignment;
using namespace std;
enum dataType {SEQUENCE,TABLE};


struct arguments
  {
    vector<unsigned> popsamples;
    unsigned nperms,outgroup;
    char *infile;
    bool haveOutgroup,bySampleSize;
    dataType type;
  };


void parseargs(int argc,char *argv[],arguments *params);
void usage(void);

int main(int argc, char *argv[])
{
  arguments params;
  parseargs(argc,argv,&params);
  unsigned totsam=accumulate(params.popsamples.begin(),params.popsamples.end(),unsigned(0));

  //vector<Fasta > data;
  //PolySites *poly;
  Sequence::PolyTable *poly = NULL;
  try
    {
     // GetData(data,params.infile);

//if sequence
      if (params.type == SEQUENCE)
	{
	  vector< Sequence::Fasta > data;
	  Sequence::Alignment::GetData(data,params.infile);
		if(Gapped(data))
        RemoveTerminalGaps(data);
	  poly = new Sequence::PolySites(data);
	  //sites_compared = Sequence::Alignment::UnGappedLength(seqs);
	}
      else if (params.type == TABLE)
	{
	  poly = new Sequence::Hudson2001;
	  ifstream in(params.infile);
	  in >> *poly;
	  in.close();
	}

 
/*      if (params.haveOutgroup == true)
        {
          //remove outgroup from the data
          vector<Fasta> newdata;
          for(unsigned i = 0 ; i < data.size() ; ++i)
            {
              if (i != params.outgroup)
                newdata.push_back(data[i]);
            }
          data.assign(newdata.begin(),newdata.end());
          newdata.clear();
        }
  *///    poly = new PolySites(data);
      
      //
      
      
    }
  catch (SeqException &e)
    {
      e.print(cerr);
      exit(10);
    }
  assert(totsam == data.size());

  int npops = params.popsamples.size();

  vector<double> weights;
  if (params.bySampleSize == true)
    {
      //weight populations by sample size
      weights.resize(npops);
      for(unsigned i = 0 ; i < params.popsamples.size() ; ++i)
        weights[i] = double(params.popsamples[i])/double(totsam);
    }
  else
    {
      //weight equally
      weights.assign(npops,1./double(npops));
    }

  //get observed values;
  double obsFST = 0.;
  try
    {
      FST obs(poly,npops,&params.popsamples[0],&weights[0]);
      obsFST = obs.HBK();
    }
  catch (SeqException &e)
    {
      cerr << e << endl;
      exit(1);
    }

  //  now, do the permutations
  Numerology::UniformDeviate Urand;
  unsigned prob=0;
  for (unsigned i = 0 ; i < params.nperms ; ++i)
    {
      random_shuffle(poly->begin(),poly->end(),Urand);
      try
        {
          FST perm(poly,npops,&(params.popsamples)[0],&weights[0]);
          if (perm.HBK() >= obsFST)
            ++prob;
        }
      catch (SeqException &e)
        {
          cerr << e << endl;
          exit(1);
        }
    }

  cout <<obsFST<<"\t"<<double(prob)/double(params.nperms)<<endl;
}

void parseargs(int argc,char *argv[],arguments *params)
{
  if(argc==1)
    {
      usage();
      exit(1);
    }
  extern int optind;
  int c;
  params->nperms = 10000;
  params->infile = NULL;
  params->haveOutgroup = false;
  params->outgroup = 0;
  params->bySampleSize = false;
  params->type = SEQUENCE;
  unsigned npop=0,dummy;
  //LEOPARD PROBLEM?
  while ((c = getopt (argc, argv, "i:h:c:n:O:s")) != -1)
    {
      switch(c)
        {
//        case 'i':
 //         params->infile = optarg;
 //         break;
 	case 'i':
	  params->infile = optarg;
	  params->type = SEQUENCE;
	  break;
	case 'h':
	  params->infile = optarg;
	  params->type = TABLE;
	  break;
        case 'c':
          npop = atoi(optarg);
	  dummy=0;
          for (unsigned k = optind ; k < optind+npop ; ++k,++dummy)
            {
              params->popsamples.push_back(atoi(argv[k]));
            }
	  optind+=dummy;
          break;
        case 'n':
          params->nperms = atoi(optarg);
          break;
        case 'O':
          params->haveOutgroup = true;
          params->outgroup = atoi(optarg);
          break;
        case 's':
          params->bySampleSize = true;
          break;
        }
    }
}

void usage(void)
{}

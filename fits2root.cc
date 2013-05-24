#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "globalConstants.h"


#include "TFile.h"
#include "TNtuple.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program produces a ROOT Ntuple from a fits image. It can handle single and two\n";
    cout << "sided images. If you are processing a two sided image you should use the \'-t\' option.\n";
    cout << "In this case the left side of the image is mirrored in the x direction.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> -<r/l>\n\n";
  cout << "\nOptions:\n";
  cout << "  -t to process a two sided image.\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -s <HDU number> for processing a single HDU \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}


int computeImage(const char *inFile, const char *outFile, vector<int> &singleHdu, const bool hasTwoSides){
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  int nhdu = 0;
  
  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&infptr, inFile, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhdu, &status);
  if (status != 0) return(status);
  
  
  for(unsigned int i=0;i<singleHdu.size();++i){
    if(singleHdu[i] > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
  }
  
  if(singleHdu.size() == 0){
    for(int i=0;i<nhdu;++i){
      singleHdu.push_back(i+1);
    }
  }
  const unsigned int nUseHdu=singleHdu.size();
  if(gVerbosity) cout << "The output file will contain " << singleHdu.size() << " of " << nhdu << " hdus availables in the input files."<< endl;
  
  fitsfile *outfptr;   /* FITS file pointers defined in fitsio.h */
  fits_create_file(&outfptr, outFile, &status);
  
  
  vector< double* > vPix;
  vector< int > vFullNCol;
  vector< int > vNLines;
  vector< string > vOHDU;
  vector< long > vTotPix;
  vector< long > vCol;
  vector< long > vRow;
  
  for (unsigned int eI=0; eI<nUseHdu; ++eI)  /* Main loop through each extension */
  {
    const unsigned int n = singleHdu[eI];
    const int nHDUsToProcess = nUseHdu;
    
    /* get input image dimensions and total number of pixels in image */
    int hdutype, bitpix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    vTotPix.push_back(totpix);
    vCol.push_back(naxes[0]);
    vRow.push_back(naxes[1]);
    double bzero;
    ffgky(infptr, TDOUBLE, "BZERO", &bzero, NULL, &status);
    if (status){
      status = 0;
      bzero = 0.0;
    }
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      continue;
    }
    
    double* outArray = new double[totpix];
    vPix.push_back(outArray);
    
    for(int i=0;i<totpix;++i) outArray[i] = 0;

    int xMin=1;
    int xMax=naxes[0];
    int yMin=1;
    int yMax=naxes[1];
    
    /* Read the images as doubles, regardless of actual datatype. */
    long fpixel[2]={xMin,yMin};
    long lpixel[2]={xMax,yMax};
    long inc[2]={1,1};
    fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, outArray, &anynul, &status);
    if (status != 0) return(status);
    
    vFullNCol.push_back(naxes[0]);
    vNLines.push_back(naxes[1]);
    
    if(gVerbosity){
      showProgress(eI,3*nHDUsToProcess);
    }
  }
  /* Close the fits files */
  fits_close_file(infptr,  &status);

  //save the root file
  ostringstream varList;
  varList << "x:y";
  long totpixNtuple=0;
  long colNtuple=0;
  long rowNtuple=0;
  int nVars=0;
  for(unsigned int e=0;e<vTotPix.size();++e){
    
    if(vTotPix[e]==1) continue;
    if(hasTwoSides){
      varList << ":L_hdu" << singleHdu[e];
      varList << ":R_hdu" << singleHdu[e];
    }
    else{
      varList << ":hdu" << singleHdu[e];
    }
    if(totpixNtuple==0){
      totpixNtuple=vTotPix[e];
      colNtuple=vCol[e];
      rowNtuple=vRow[e];
    }
    
    if(colNtuple!=vCol[e] || rowNtuple!=vRow[e]){
      cerr << "\n\nError, the images don't have the same size!\n\n";
      cout << vTotPix[e] << "\t" << colNtuple << "\t" << rowNtuple << "\t" << vCol[e] << "\t" << vRow[e] << "\n";
      return -1000;
    }
    ++nVars;
  }
  
  if(hasTwoSides){
    colNtuple/=2;
    totpixNtuple/=2;
  }
  
  
  TFile outF( outFile, "RECREATE");
  TNtuple nt("data","data",varList.str().c_str());
  
  Float_t *varArr;
  if(hasTwoSides){
    varArr = new Float_t[2+2*nVars];
  }
  else{
    varArr = new Float_t[2+nVars];
  }
  
  const unsigned int nExt=vPix.size();
  for(long i=0;i<totpixNtuple;++i){
    
    long x = i%colNtuple;
    long y = i/colNtuple;
    
    varArr[0] = x;
    varArr[1] = y;
    
    for(unsigned int e=0;e<nExt;++e){
      if(vTotPix[e]==0) continue;
      if(hasTwoSides){
        varArr[2+2*e]   = vPix[e][colNtuple - x - 1 + 2*y*colNtuple];
        varArr[2+2*e+1] = vPix[e][colNtuple + x + 2*y*colNtuple];
      }
      else{
        varArr[2+e] = vPix[e][i];
      }
    }
    nt.Fill(varArr);
    
    if(gVerbosity){
      if(y%500==0) showProgress((totpixNtuple+i*2)/3,totpixNtuple);
    }
  }
  
  nt.Write();
  outF.Close();
  
  
  /* clean up */
  for(unsigned int i=0;i<vPix.size();++i){
    delete[] vPix[i];
  }
  
  if(gVerbosity){
    showProgress(1,1);
  }
  
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], vector<int> &singleHdu, string &inFile, string &outFile, bool &hasTwoSides){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  singleHdu.clear();
  int opt=0;
  while ( (opt = getopt(argc, argv, "i:o:s:qthH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      singleHdu.push_back(atoi(optarg));
      break;
    case 't':
      hasTwoSides = true;
      break;
    case 'q':
      gVerbosity = 0;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }
  
  inFile="";
  
  if(argc-optind==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  else if(argc-optind>1){
    cerr << red << "Error: more than one input file provided!\n\n" << normal;
    return 1;
  }
  
  inFile=argv[optind];
  if(!fileExist(inFile.c_str())){
    cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string outFile;
  string inFile;
  vector<int> singleHdu;
  bool hasTwoSides=false;
  
  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFile, outFile, hasTwoSides);
  if(returnCode!=0){ 
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist. " << normal;
    deleteFile(outFile.c_str());
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following file:\n" << normal;
    cout << "\t" << inFile << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  int status = computeImage( inFile.c_str(),  outFile.c_str(), singleHdu, hasTwoSides);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}

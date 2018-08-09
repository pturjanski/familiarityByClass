#include "smrUnn.h"


/**
 * Finds maximum between two numbers.
 */
int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}
 
 
 
/**
 * Finds minimum between two numbers.
 */
int min(int num1, int num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}


smrUnn::smrUnn(const string& inputFastaFolderName){
    bool ok;
    vector<FastaElement> fs;
    ok=Tools::loadFastaFilesFromFolder(inputFastaFolderName, fs);
    if(!ok) cout << "Class mr->loadFastaFilesFromFolder: Error!" << endl;
    
    //1. total=#aminoacids
    unsigned int total=0;
    for(unsigned int iProtein=0;iProtein<fs.size();iProtein++){
        total+=fs[iProtein].seq.size();
    }
    //2. memory assignment = total+#separatorSymbols+2  (2=endTerminatorSymbol$ and endTerminatorSymbol|)
    n=total+fs.size()+2;
    T = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    //3. Copy aminoacids+separatorSymbols to T
    unsigned int iPos=0;
    for(unsigned int iProtein=0;iProtein<fs.size();iProtein++){
        for(unsigned int iSeq=0;iSeq<fs[iProtein].seq.size();iSeq++){
            T[iPos]=fs[iProtein].seq[iSeq];
            iPos++;
        }
        T[iPos]='+';
        iPos++;
   }
   T[iPos]='$';
   iPos++;
   T[iPos]='|'; // > Alphabet (used by l_intervals algorithm)

   computeSA();
   computeLCP();
   computePatterns();
}
    
/* Construct the suffix array. */
void smrUnn::computeSA(){
    SA  = (int *)malloc((size_t)n * sizeof(int));
    if(SA == NULL) {
        cout << "Class mr->computeSA: Cannot allocate memory." << endl;
        exit(EXIT_FAILURE);
    }
    if(sais(T, SA, (int)n) != 0) {
        cout << "Class mr->computeSA->sais: Cannot allocate memory." << endl;
        exit(EXIT_FAILURE);
    }
}

void smrUnn::computeLCP(){
    // require SA computed
    // ... naive LCP:
    LCP = (int *)malloc((size_t)n * sizeof(int));

    if(LCP == NULL) {
        cout << "Class mr->computeLCP: Cannot allocate memory." << endl;
        exit(EXIT_FAILURE);
    }

    for (int iSA = 1; iSA < n; ++iSA) {
        int l = 0;
        while ((T[SA[iSA]+l]==T[SA[iSA-1]+l]) &&  T[SA[iSA]+l]!='+') ++l;
        LCP[iSA] = l;
    }

    // In order to avoid border cases in repeat classification, assigns
    //    ... LCP[0]=0
    LCP[0]=0;
}   



smrUnn::~smrUnn(){
  free(T);
  free(LCP);
  free(SA);
}


void smrUnn::computePatterns(){
  //require: SA, LCP

  // allocate memory 
  MaxRep  = (int *)malloc((size_t)n * sizeof(int));
  if(MaxRep == NULL) {
        fprintf(stderr, "Class smrUnn->computeSMR: Cannot allocate memory.\n");
        exit(EXIT_FAILURE);
  }

  // Result vector initialization
  for(unsigned int idx=0;idx<(unsigned int)n;idx++) MaxRep[idx]=n;
  // SMR and NN finding
  for(unsigned int idx=0;idx<(unsigned int)n;idx++) {
    int lcpValue=max(LCP[idx],LCP[idx+1]);
    if(lcpValue>0){ // (Only repetition analysis)
        int patternEndPosition   = SA[idx]+lcpValue-1;
        int patternStartPosition = min(MaxRep[patternEndPosition],SA[idx]);
        MaxRep[patternEndPosition]=patternStartPosition;
    }
  }
  
  // Insert patterns into SMR_NN set
  for(unsigned int idx=0;idx<n;idx++){
    if(MaxRep[idx]!=n){
        string pattern="";
        for(unsigned int ipos=MaxRep[idx];ipos<=idx;ipos++) pattern+=T[ipos];
        SMR_NN.insert(pattern);
    }
  }
  
}



#include "mr.h"

/**
 * sets union function
 */
template<typename S>
S union_sets(const S& s1, const S& s2)
{
     S result = s1;

     result.insert(s2.begin(), s2.end());

     return result;
}

mr::mr(const string& inputFastaFolderName){
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
void mr::computeSA(){
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

void mr::computeLCP(){
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



mr::~mr(){
  free(T);
  free(LCP);
  free(SA);
}


void mr::addOccurrence(map<char, unsigned int> &map1, char e)
{
    if(map1.find(e) != map1.end()){
        // (+) is a distinguished symbol. + != + . At most it can only have one instance.
        if(e!='+') map1[e]++;  
    }    
    else{
        map1.insert(map<char,unsigned int>::value_type(e,1));
    }
}


void mr::l_intervals(vector<unsigned int> &vector_l, vector<unsigned int> &vector_i, vector<unsigned int>& vector_j){
    
    vector<unsigned int> interval_l;
    vector<unsigned int> interval_i;
    vector<unsigned int> interval_j;
    
    interval_l.push_back(0);
    interval_i.push_back(0);
    interval_j.push_back(0); // NULL
    for(unsigned int i=1;i<=n;i++){
        unsigned int lb = i - 1;
        while (LCP[i] < (int)interval_l[interval_l.size()-1]){
            interval_j[interval_j.size()-1]=i-1;
            
            vector_l.push_back(interval_l[interval_l.size()-1]);
            vector_i.push_back(interval_i[interval_i.size()-1]);
            vector_j.push_back(interval_j[interval_j.size()-1]);
            
            
            lb = interval_i[interval_i.size()-1];
            interval_l.pop_back();
            interval_i.pop_back();
            interval_j.pop_back();
            
        }
        if(LCP[i] > (int) interval_l[interval_l.size()-1]){
            interval_l.push_back(LCP[i]);
            interval_i.push_back(lb);
            interval_j.push_back(0); // NULL
        }
    }
}




void mr::computePatterns(){
    // require SA, LCP
    
    // 1. Compute Intervals
    vector<unsigned int> ls;
    vector<unsigned int> is;
    vector<unsigned int> js;
    l_intervals(ls, is, js);
        
    // 2. Parameter set
    LCP[0]=0; // Ojo, es de prueba!    

    // 3. each pattern classification (l-interval patterns)
    for(unsigned int it=0;it<is.size();it++){ 
        unsigned int l=ls[it];
        unsigned int i=is[it];
        unsigned int j=js[it];
        
        // pattern
        string pattern="";
        for(unsigned int isuff=0;isuff<l;isuff++) pattern+=(T[SA[i]+isuff]); 
        // variables initialization
        map<char, unsigned int> LCTT, RCTT; //empty dictionaries
        bool is_MR=false;
        char sigma = T[SA[i]-1];
        // left & right context count
        for(unsigned int k=i;k<=j;k++){
            char sigma_l=T[SA[k]-1]; // left  context
            char sigma_r=T[SA[k]+l]; // right context
            addOccurrence(LCTT, sigma_l);   
            addOccurrence(RCTT, sigma_r);
	    // l-interval always contains at least one element 
            // that is not-right extensible (by definition).
            // We only need to check if it's not left-extensible            
            if(sigma!=sigma_l || sigma=='+') is_MR = true;
        }

        // is it SMR, NN or NE?
        if(is_MR){    
            // MR is a SMR?
            bool allLessOrEqOne = true;
            for ( std::map<char,unsigned int>::iterator it = LCTT.begin(); it != LCTT.end() && allLessOrEqOne; ++it){
                allLessOrEqOne = (allLessOrEqOne && (it->second<=1));
            }
            for ( std::map<char,unsigned int>::iterator it = RCTT.begin(); it != RCTT.end() && allLessOrEqOne; ++it){
                allLessOrEqOne = (allLessOrEqOne && (it->second<=1));
            }
            
            if(allLessOrEqOne){
                // There is not pair of occurrences that has the same left or right-context
                // Then, they can't be extended
                SMR.insert(pattern); // It's SMR
            }
            else{
                // is it  NN or NE?
                // Find an occurrence that is not right and left extensible
                // If it exist then it is a NN 
                bool isNested = true;
                for(unsigned int k=i;(k<=j) && isNested;k++){ // for each occurrence ...
                    char sigma_l=T[SA[k]-1]; // left  context
                    char sigma_r=T[SA[k]+l]; // right context
                    if( (LCTT[sigma_l]<=1) && (RCTT[sigma_r]<=1) ){
                        // This occurrence is not right and left extensible
                        isNested = false;
                    }
                }
                if(isNested) NE.insert(pattern); // It is a NE
                                                    // (all pattern occurrences are nested)
                else         NN.insert(pattern); // It is a NN
                                                    // (at least one pattern occurrence is non-nested)
            }
        }
    }
}


set<string> mr::SMR_NN_NE() const{
 set<string> SMR_NN=union_sets(SMR, NN);
 return union_sets(SMR_NN, NE);
}

set<string> mr::SMR_NN() const{
 return union_sets(SMR, NN);
}

set<string> mr::SMR_NE() const{
 return union_sets(SMR, NE);
}

set<string> mr::NN_NE() const{
 return union_sets(NN, NE);
}
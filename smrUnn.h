#ifndef SMRUNN_H
#define SMRUNN_H
#include <string>
#include <set>
#include "sais.h"
#include "Tools.h"
#include "FastaElement.h"



using namespace std;


class smrUnn {

    private:
        unsigned char *T;
        int *SA;
        int *LCP;  
        int *MaxRep;
        long n;

        // Methods
        void computeSA();
        void computeLCP();

        void computePatterns();
        
    public:
        // Patterns
        set<string> SMR_NN;


        // Constructor/ Destructor
        smrUnn(const string& inputFastaFolderName);
        ~smrUnn();
        
};

#endif // SMRUNN_H



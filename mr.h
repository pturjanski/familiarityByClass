#ifndef MR_H
#define MR_H
#include <string>
#include <set>
#include <map>
#include <vector>

#include "sais.h"
#include "Tools.h"
#include "FastaElement.h"



using namespace std;


class mr {

    private:
        unsigned char *T;
        int *SA;
        int *LCP;  
        long n;

        // Methods
        void repeatClassification(const unsigned int l, const unsigned int i, const unsigned int j) const;
        void l_intervals(vector<unsigned int> &l, vector<unsigned int> &i, vector<unsigned int>& j);

        void computeSA();
        void computeLCP();

        void computePatterns();
        
        void addOccurrence(map<char, unsigned int> &map1, char e);
    public:
        // Patterns
        set<string> NE;
        set<string> NN;
        set<string> SMR;

        set<string> SMR_NN_NE() const;
        set<string> SMR_NN() const;
        set<string> SMR_NE() const;
        set<string> NN_NE() const;

        // Constructor/ Destructor
        mr(const string& inputFastaFolderName);
        ~mr();
        
};

#endif // MR_H



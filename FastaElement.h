#ifndef FASTAELEMENT_H
#define FASTAELEMENT_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "Trie.h"



using namespace std;

class FastaElement {
    private:
        string getFilename(const string s) const;
   
    
    public:

	string name;
	string description;
	string seq; 
	string family;

        bool loadFromFastaFile(const string filename);
        double familiarity_10(const Trie*  trie) const;
        double coverage(const Trie*  trie, const unsigned int minimunRepeatLength) const;
        
        bool operator==(const FastaElement &f2) const;  

	// Constructor
	FastaElement();
        
};

#endif // FASTAELEMENT_H

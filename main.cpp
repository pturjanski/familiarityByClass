#include <iostream>
#include "Trie.h"

using namespace std;

// Test program
int main()
{
    Trie* trie = new Trie();
    trie->addWord("");
    trie->addWord("Hello");
    trie->addWord("Balloon");
    trie->addWord("Ball");

    
    vector<string> success;
    success=trie->prefixMatching("ABalloonPie");
  
    for(unsigned int i=0;i<success.size();i++){
        cout << success[i] << endl;
    }
    
    
//     if ( trie->searchWord("Hell") )
//         cout << "Found Hell" << endl;
// 
//     if ( trie->searchWord("Hello") )
//         cout << "Found Hello" << endl;
// 
//     if ( trie->searchWord("Helloo") )
//         cout << "Found Helloo" << endl;
// 
//     if ( trie->searchWord("Ball") )
//         cout << "Found Ball" << endl;
// 
//     if ( trie->searchWord("Balloon") )
//         cout << "Found Balloon" << endl;

    delete trie;
}

#ifndef TRIE_H
#define TRIE_H

#include <vector>
#include <set>
#include <string>
#include "Node.h"
using namespace std;

class Trie {
public:
    Trie();
    ~Trie();

    void addWord(const string& s);
    void addSetOfWords(const set<string>& s);

    bool searchWord(const string& s) const;
    vector<string> prefixMatching(const string& t) const;
    vector<string> prefixMatchingFromStartPosition(const string& t, const unsigned int p) const;
    unsigned int LongestPatternMatchingFromStartPosition(const string& t, const unsigned int p) const;

private:
    void deleteTrie(Node* r);
    Node* root;
};
#endif // TRIE_H

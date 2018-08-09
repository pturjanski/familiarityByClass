#include "Trie.h"

Trie::Trie(){
    root=new Node();
}

Trie::~Trie(){
    // Free memory
    deleteTrie(root);
}

void Trie::deleteTrie(Node* r){
    if(r==NULL) return;
    for(unsigned int i=0; i<r->children().size(); i++)
        deleteTrie(r->children()[i]);
    delete r;
    r = NULL;
}

void Trie::addWord(const string& s){
    Node* current = root;

    if(s.length()==0){
        current->setWordMarker(); // an empty word
        return;
    }

    for (unsigned int i=0;i<s.length();i++){        
        Node* child=current->findChild(s[i]);
        if(child!=NULL){
            current=child;
        }
        else{
            Node* tmp=new Node();
            tmp->setContent(s[i]);
            current->appendChild(tmp);
            current=tmp;
        }
        if(i==s.length()-1) current->setWordMarker();
    }
}


void Trie::addSetOfWords(const set<string>& s){
    for(set<string>::iterator si = s.begin(); si!=s.end(); si++) addWord(*si);
}


bool Trie::searchWord(const string& s) const{
    Node* current=root;

    while(current!=NULL){
        for (unsigned int i=0; i<s.length();i++){
            Node* tmp=current->findChild(s[i]);
            if(tmp==NULL) return false;
            current=tmp;
        }

        if (current->wordMarker()) return true;
        else return false;
    }
    return false;
}



vector<string> Trie::prefixMatching(const string& t) const{
    
    return prefixMatchingFromStartPosition(t,0);
}

vector<string> Trie::prefixMatchingFromStartPosition(const string& t, const unsigned int p) const{
    
    vector<string> res;
    Node* current = root;
    // check for the empty pattern
    if(current->wordMarker()) res.push_back("");

    // case t is a valid position
    if (t.size()>p){ 
        // remaining patterns
        unsigned int i=p;
        while(current!=NULL && i<t.length()){
            Node* tmp=current->findChild(t[i]);
            current=tmp;
            if (tmp!=NULL && current->wordMarker()) res.push_back(t.substr(p,i-p+1)); // pattern
            i++;
        }
    }
    return res;
}

unsigned int Trie::LongestPatternMatchingFromStartPosition(const string& t, const unsigned int p) const{
    
    unsigned int res=0;
    Node* current = root;

    // case t is a valid position
    if (t.size()>p){ 
        // remaining patterns
        unsigned int i=p;
        while(current!=NULL && i<t.length()){
            Node* tmp=current->findChild(t[i]);
            current=tmp;
            if (tmp!=NULL && current->wordMarker()) res=i-p+1; // pattern size
            i++;
        }
    }
    return res;
}

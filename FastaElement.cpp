#include "FastaElement.h"

FastaElement::FastaElement(){
	name="";
	description="";
	seq="";
	family="";
}


bool FastaElement::operator==(const FastaElement&f2) const{
	return (seq==f2.seq);
};




double FastaElement::coverage(const Trie*  trie, const unsigned int minimunRepeatLength) const{
    // vector initialization
    vector<bool> affectedPositions;
    for(unsigned int i=0;i<seq.size();i++) affectedPositions.push_back(false);

    // mark with true affected positions
    for(unsigned int i=0;i<seq.size();i++){
        unsigned int maxPattern = trie->LongestPatternMatchingFromStartPosition(seq,i);
        if(maxPattern>=minimunRepeatLength){
            for(unsigned int k=i;k<(i+maxPattern);k++) affectedPositions[k]=true;
        }
    }
    
    //count affected positions
    double ret=0;
    for(unsigned int i=0;i<affectedPositions.size();i++){
        //fprintf(stderr," %d \n",affectedPositions[i]?1:0);
        if(affectedPositions[i]) ret=ret+1.0;
    }
    return (ret/seq.size());
}
  

double FastaElement::familiarity_10(const Trie*  trie) const{
    double sum=0;
    for(unsigned int i=1;i<10;i++) sum+=coverage(trie, i);
    return ( (1+coverage(trie, 10))/2.0 + sum );
}




bool FastaElement::loadFromFastaFile(const string filename){

  vector<FastaElement> fs;
  // File Open 
  ifstream is(filename.c_str());

  if(!is.good()){
      cout << "Error! function " << __FUNCTION__ << " cannot open " << filename << endl;
      return false;
  }
  else{
    
      // Load protein (if there are more than one, all of them are loaded into fs)
      char c;
      if(!is.eof()) {
	    is >> c;
	    if(c!='>') {
		    cout << "Error! function " << __FUNCTION__ << " ! Multifasta file has an incorrect format." << endl;
		    return false;
	    }
      }
      while(!is.eof()){
	    FastaElement f;
	    f.name=getFilename(filename); // delete path and extension from filename
	    getline(is,f.description);
	    is >> c;
	    while((!is.eof()) && (c!='>')){
		    f.seq+=c;
		    is >> c;
	    }
	    f.family="";
	    fs.push_back(f);
      }
      is.close();
  }
  
  
  // Only one protein?
  if(fs.size()!=1){
      cout << "Error! Function " << __FUNCTION__ << " . " << filename << " is not fasta file. It's a MULTIFASTA! This file has " << fs.size() << " proteins. We allow only one sequence." << endl;
      return false;
  }
    
  // f is the unique protein in the file
  FastaElement f=fs[0];
    
  // assign f to this
  name=f.name;
  description=f.description;
  seq=f.seq;
  family="";  
  
  return true;

}

string FastaElement::getFilename(const string s) const{
  ///////////////////////////////////////////////////  
  //     delete path and extension from filename   // 
  ///////////////////////////////////////////////////
  string fileName;
  // fileName bars counter
  int barsCounter=0;
  for(int icb=0;icb<(int)s.size();icb++) if(s[icb]=='/') barsCounter++;
  istringstream ss(s);
  // Remove all before bars
  for(int icb=0;icb<barsCounter;icb++) getline(ss,fileName,'/');
  // Remove extension
  getline(ss,fileName,'.');
  return fileName;
}
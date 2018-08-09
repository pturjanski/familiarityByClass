#include "Tools.h"



string Tools::getFilename(const string s){
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




bool Tools::loadMultifastaFile(const string filename, vector<FastaElement>& fs){
  bool ret=true;
  // File open
  ifstream is(filename.c_str());
  if(!is.good()){
      cout << "Error! function " << __FUNCTION__ << " cannot open " << filename << endl;
      ret=false;
  }
  else{
    
      // load protein 
      char c;
      if(!is.eof()) {
	    is >> c;
	    if(c!='>') {
		    cout << "Error! function " << __FUNCTION__ << " ! Multifasta file has an incorrect format." << endl;
		    return 1;
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
  return ret;

}



bool Tools::loadFastaFilesFromFolder(const string folder, vector<FastaElement>& fs){

  path p(folder); 

  try
  {
    if (exists(p))    // does p actually exist?
    {
      if (is_regular_file(p)){        // is p a regular file?
        cout << "Error! function " << __FUNCTION__ << " !" << folder << " not is a directory!" << endl;
        return false;
      }
      else if (is_directory(p))      // is p a directory?
      {
        
        // for each file in folder
        for (directory_iterator itr(p); itr!=directory_iterator(); ++itr)
        {
            string proteinFileName = itr->path().filename().string();
            string proteinFileNameWithPath = itr->path().string();
            if( (proteinFileName.size()>=6) && (proteinFileName.substr(proteinFileName.size()-6,6)==".fasta") ) // is fasta file?
            {
                vector<FastaElement> aux_ps;
                bool ok;
                ok=Tools::loadMultifastaFile(proteinFileNameWithPath, aux_ps);
                if(!ok) return false;
                if(aux_ps.size()!=1){
                    cout << "Error! Function " << __FUNCTION__ << " . " << proteinFileNameWithPath << " is not fasta file. It's a MULTIFASTA! This file has " << aux_ps.size() << " proteins. The problem resides at the moment to assign the name (it always assigns the filename)." << endl;
                    return false;
                }
                else{
                    fs.push_back(aux_ps[0]);
                }
            }
            else{
                //not fasta file
                cout << "NOT FASTA: " << itr->path() << endl;
            }
        }
      }
      else{
        cout << "Error! function " << __FUNCTION__ << " !" << folder << " exists, but is neither a regular file nor a directory\n";
        return false;
      }
    }
    else{
      cout << "Error! function " << __FUNCTION__ << " !" << folder << " does not exist!" << endl;
      return false;
    }
  }

  catch (const filesystem_error& ex)
  {
    cout << "Error! function " << __FUNCTION__ << " !" << endl;
    cout << ex.what() << endl;
    return false;
  }
    
  
  return true;

}








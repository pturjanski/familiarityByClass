#include <iostream>
#include <iterator>
#include  "mr.h"
#include  "smrUnn.h"
#include  "Trie.h"


void familiarityAndCoverageDetails(const FastaElement& sequence, const set<string>& setOfpatterns, double& familiarity, vector<double>& coverageDetails){
    // Trie from patterns (allows search optimization)
    Trie* trie = new Trie();
    trie->addSetOfWords(setOfpatterns);
    // compute familiarity        
    familiarity = sequence.familiarity_10(trie); 
    // compute coverage details
    for(unsigned int minimunPatternLength=0;minimunPatternLength<=10;minimunPatternLength++) coverageDetails.push_back(sequence.coverage(trie, minimunPatternLength));
    // memory free
    delete trie;
    
}

void computeAndSave(const FastaElement& sequence, const set<string>& setOfpatterns, const string method, ostream& osFamiliarity, ostream& osCoverage){
    double familiarity;
    vector<double> coverageDetails;
    // Compute familiarity and coverage details
    familiarityAndCoverageDetails(sequence, setOfpatterns, familiarity, coverageDetails);
    // save familiarity
    osFamiliarity << method;
    osFamiliarity << "," << familiarity << endl;
    // save coverage details
    osCoverage << method;
    for(unsigned int minimunPatternLength=0;minimunPatternLength<=10;minimunPatternLength++) osCoverage << "," << coverageDetails[minimunPatternLength];
    osCoverage << endl;
    
}

bool saveSet(const string outputFilename, const set<string>& s){
    // ... Open output file
    std::ofstream osData(outputFilename.c_str());
    if(!osData.good()){
        cout << "Error! function " << __FUNCTION__ << " cannot open " << outputFilename << endl;
        return false;
    }
    // ... File header 
    osData << "Pattern" << endl;
    // ... Save patterns
    for (std::set<string>::iterator it = s.begin(); it != s.end(); ++it){
        osData << *it << endl; 
    }
    // ... Close file
    osData.close();
    
    return true;
}


int main(int argc, const char *argv[]) {

  // paramaeter checking
    if(argc!=4){
                cout<<"usage example: " << argv[0] << " sequence.fasta familyDatasetFolder outputFolder" << endl;
                return 0;
    }


    bool ok;

    // 1. Parameters
    // =============
    //a. sequenceFileName
    string sequenceFileName =  argv[1];

    //b. familyDatasetFolder
    string familyDatasetFolder =  argv[2];


    //c. outputFolder
    string outputFolder = argv[3];

    cout << endl << endl << endl; 
    cout << "Parameters" << endl;
    cout << "==========" << endl;
    cout << " 1. sequence            : " << sequenceFileName    << endl;
    cout << " 2. familyDatasetFolder : " << familyDatasetFolder << endl;
    cout << " 3. outputFolder        : " << outputFolder        << endl;
    
    cout << endl << endl;
    
    // 1. Load sequence
    // ================
    FastaElement sequence;
    ok=sequence.loadFromFastaFile(sequenceFileName);
    if(!ok) return 1;
    
    // 2. Compute Patterns
    // ===================
    //////////////////////////////////////////////////////////////////////////////////
    ////////                                                                  ////////
    //////// Using Listing 1 "MR Identification and Classification Algorithm" ////////
    ////////                                                                  ////////
    //////////////////////////////////////////////////////////////////////////////////
    mr mrSet(familyDatasetFolder); 
    
    // 3. Save sets
    // ============
    ok=saveSet(outputFolder+"SMR_NN_NE.csv", mrSet.SMR_NN_NE());
    if(!ok) return 1;
    ok=saveSet(outputFolder+"SMR.csv"      , mrSet.SMR);
    if(!ok) return 1;
    ok=saveSet(outputFolder+"NN.csv"       , mrSet.NN );
    if(!ok) return 1;
    ok=saveSet(outputFolder+"NE.csv"       , mrSet.NE );
    if(!ok) return 1;
    ok=saveSet(outputFolder+"SMR_NN.csv"   , mrSet.SMR_NN());
    if(!ok) return 1;
    ok=saveSet(outputFolder+"SMR_NE.csv"   , mrSet.SMR_NE());
    if(!ok) return 1;
    ok=saveSet(outputFolder+"NN_NE.csv"    , mrSet.NN_NE() );
    if(!ok) return 1;

    // 4. Open output files
    // ====================
    //// a. familiarity
    string outputFamiliarityFilename=outputFolder+"familiarity.csv";
    std::ofstream osFamiliarity(outputFamiliarityFilename.c_str());
    if(!osFamiliarity.good()){
        cout << "Error! File " << outputFamiliarityFilename << " cannot be opened " << endl;
        return 1;
    }
    //// ... header of file (sequence name)
    osFamiliarity << "Filename=" << sequence.name << endl;

    //// b. coverage
    string outputCoverageFilename=outputFolder+"coverage.csv";
    std::ofstream osCoverage(outputCoverageFilename.c_str());
    if(!osCoverage.good()){
        cout << "Error! File " << outputCoverageFilename << " cannot be opened " << endl;
        return 1;
    }
    //// ... header of file
    osCoverage << "Method";
    for(unsigned int minimunPatternLength=0;minimunPatternLength<=10;minimunPatternLength++) osCoverage << ",coverage(seq,minimunPatternLength(patterns," << minimunPatternLength << "))";
    osCoverage << endl;

    // 5. Compute Familiarity and Details
    // ==================================
    
    // Familiarity and coverage details with MR (SMR U NN U NE)
    computeAndSave(sequence,mrSet.SMR_NN_NE(),"SMR_NN_NE (MR)",osFamiliarity, osCoverage);
    // Familiarity and coverage details with SMR
    computeAndSave(sequence,mrSet.SMR        ,"SMR"           ,osFamiliarity, osCoverage);
    // Familiarity and coverage details with NN
    computeAndSave(sequence,mrSet.NN         ,"NN"            ,osFamiliarity, osCoverage);
    // Familiarity and coverage details with NE
    computeAndSave(sequence,mrSet.NE         ,"NE"            ,osFamiliarity, osCoverage);
    // Familiarity and coverage details with SMR U NN
    computeAndSave(sequence,mrSet.SMR_NN()   ,"SMR_NN"        ,osFamiliarity, osCoverage);
    // Familiarity and coverage details with SMR U NE
    computeAndSave(sequence,mrSet.SMR_NE()   ,"SMR_NE"        ,osFamiliarity, osCoverage);
    // Familiarity and coverage details with NN U NE
    computeAndSave(sequence,mrSet.NN_NE()    ,"NN_NE"         ,osFamiliarity, osCoverage);
    
    
    
    // 6. Close Files
    // ==============
    osFamiliarity.close();
    osCoverage.close();


    // 7. Compute Patterns
    // ===================
    //////////////////////////////////////////////////////////////////////////////////
    ////////                                                                  ////////
    //////// Using Listing 3 "Efficient SMR and NN Identification Algorithm"  ////////
    ////////                                                                  ////////
    //////////////////////////////////////////////////////////////////////////////////
    smrUnn smrUnnSet(familyDatasetFolder); 

    
    // 8. Save Patterns files
    // ======================
    ok = saveSet(outputFolder+"SMR_NN_Listing3.csv", smrUnnSet.SMR_NN);  
    if(!ok) return 1;
 
    return 0;

  
}



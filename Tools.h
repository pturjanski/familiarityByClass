#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>

#include "FastaElement.h"

using namespace std;
using namespace boost::filesystem;


class Tools{
    private:
        // fasta file loader
        static bool loadMultifastaFile(const string filename, vector<FastaElement>& fs);
        // data manipulation
        static string getFilename(const string s);

    public:
        // file loader
        static bool loadFastaFilesFromFolder(const string folder, vector<FastaElement>& fs);

};
#endif // TOOLS_H
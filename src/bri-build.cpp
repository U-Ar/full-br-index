#include <chrono>
#include <iostream> 
#include <string>

#include "br_index_full.hpp"
#include "utils.hpp"

using namespace std;
using namespace bri;

string out_basename = string();
string input_file = string();
bool sais = true;
long length = 8;

void help(){
	cout << "bri-build: builds the fully functional bidirectional r-index." << endl;
    cout << "           extension .brif is automatically added to output index file." << endl << endl;
	cout << "Usage: bri-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use <basename> as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
    cout << "   -divsufsort          use divsufsort algorithm to build the BWT (fast, 7.5n Bytes of RAM). By default,"<<endl;
	cout << "                        SE-SAIS is used (about 4 times slower than divsufsort, 4n Bytes of RAM)."<<endl;
    cout << "   -bl <length>         parameter (8 by default)" << endl;
    cout << "                        For contractions of the patterns shorter than <length>, " << endl;
    cout << "                        shortcut is applied." << endl;
    cout << "                        The contraction speed improves with higher bl in exchange of the index size." << endl;
    cout << "   <input_file_name>    input text file." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-o") == 0) 
    {

		if(ptr >= argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}
    else if (s.compare("-divsufsort") == 0 )
    {
		sais = false;
	}
    else if (s.compare("-bl") == 0)
    {
        if(ptr >= argc-1){
			cout << "Error: missing parameter after -bl option." << endl;
			help();
		}

		char* e;
        length = strtol(argv[ptr],&e,10);

        if(*e != '\0' || length <= 0){
            cout << "Error: illegal or nonpositive value after -bl option." << endl;
            help();
        }

		ptr++;
    }
    else
    {
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv) 
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

    out_basename = string();
    input_file = string();

    int ptr = 1;

    if (argc < 2) help();

    while (ptr < argc-1) parse_args(argv, argc, ptr);

    input_file = string(argv[ptr]);

    if (out_basename.compare("") == 0)
        out_basename = string(input_file);
    
    string idx_file = out_basename;

    cout << "Building fully functional br-index with bl=" << length << " on " << input_file << endl; 
    cout << "Index will be saved to " << idx_file << ".brif" << endl;

    string input;

    {
        std::ifstream fs(input_file);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
    }

    br_index_full<> idx(input,(unsigned long)length,sais);
    idx.save_to_file(idx_file);
    
    auto t2 = high_resolution_clock::now();

    ulint total = duration_cast<duration<double, std::ratio<1>>>(t2-t1).count();
    cout << "Build time: " << get_time(total) << endl;

}
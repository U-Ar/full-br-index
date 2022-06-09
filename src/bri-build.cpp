#include <chrono>
#include <iostream> 
#include <string>

#include "br_index.hpp"
#include "br_index_limited.hpp"
#include "utils.hpp"

using namespace std;
using namespace bri;

string out_basename = string();
string input_file = string();
bool sais = true;
bool limited = false;

void help(){
	cout << "bri-build: builds the full bidirectional r-index. Extension .bri is automatically added to output index file" << endl << endl;
	cout << "Usage: bri-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use <basename> as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
    cout << "   -limited             build limited version of br-index optimized for finding MEMs." << endl;
    cout << "                        space-efficient, but left-extension and right-contraction are not supported." << endl;
	cout << "   -divsufsort          use divsufsort algorithm to build the BWT (fast, 7.5n Bytes of RAM). By default,"<<endl;
	cout << "                        SE-SAIS is used (about 4 time slower than divsufsort, 4n Bytes of RAM)."<<endl;
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
    else if (s.compare("-limited") == 0)
    {
        limited = true;
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

    cout << "Building br-index of input file " << input_file << endl;
    if (!limited)
    {
        cout << "Index will be saved to " << idx_file << ".bri" << endl;
    }
    else {
        cout << "Index will be saved to " << idx_file << ".bril" << endl;
    }

    string input;

    {
        std::ifstream fs(input_file);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
    }


    if (!limited)
    {
        br_index<> idx(input,sais);
        idx.save_to_file(idx_file);
    }
    else 
    {
        br_index_limited<> idx(input,sais);
        idx.save_to_file(idx_file);
    }
    
    auto t2 = high_resolution_clock::now();

    ulint total = duration_cast<duration<double, std::ratio<1>>>(t2-t1).count();
    cout << "Build time: " << get_time(total) << endl;

    out.close();
}
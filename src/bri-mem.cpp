#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index_full.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check = string();
long length = 0;

void help()
{
	cout << "bri-mem: locate maximal exact matches of the input patterns. " << endl;

	cout << "Usage: bri-mem [options] <index> <pattern>" << endl;
	cout << "   (not yet implemented) -c <text>    check correctness of each pattern occurrence on this text file" << endl;
	cout << "   -bl <length>         parameter (index file's specification by default)" << endl;
    cout << "                        For contractions of the patterns shorter than <length>, " << endl;
    cout << "                        shortcut is applied." << endl;cout << "   <index>      index file (with extension .brif)" << endl;
	cout << "   <pattern>    text file of the pattern you want to search for substrings" << endl;

	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-c") == 0) 
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;
    
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

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

template<class T>
void locate_all(ifstream& in, string pattern_file)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    string text;
    bool c = false;

    if (check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    cout << "Loading fully functional br-index ... " << flush;

    auto t1 = high_resolution_clock::now();

    T idx;

    if (length) idx.load(in,(int)length);
    else idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "done." << endl;

    cout << "searching maximal exact matches ..." << endl;
    
    
    

    ifstream ifs(pattern_file);

    stringstream ssp;
    ssp << ifs.rdbuf();
    string p = ssp.str();


    ulint m = p.size();

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    ulint tot_time = 0;



    t3 = high_resolution_clock::now();

    idx.maximal_exact_match(p);

    t4 = high_resolution_clock::now();

    tot_time += duration_cast<microseconds>(t4-t3).count();


    /*if (c) // check occurrences
    {
        cout << "length of MEM : " << res.first << endl;
        cout << "number of MEM occs : " << occs.size() << endl;
        ulint maxlen = res.first;
        for (auto pair: res.second)
        {
            ulint offset = pair.first;
            auto occs = idx.locate_sample(pair.second);
            for (auto o: occs)
            {
                for (ulint i = 0; i < maxlen; ++i)
                {
                    if (text[o+i] != p[i+offset])
                    {
                        cout << "wrong occurrence  occ: " << o << " offset on P: " << offset << endl;
                    }
                }
            }
        }
    }*/

    

    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    cout << "Pattern length m = " << m << endl;
    
    // cout << "LF-mapping time: " << count_time << " microseconds" << endl;
    // cout << "Phi        time: " << locate_time << " microseconds" << endl;
    cout << "Computation time : " << tot_time << " microseconds" << endl;
    
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc-2) parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string patt_file(argv[ptr+1]);

    char* e;

    ifstream in(idx_file);


    locate_all<br_index_full<> >(in, patt_file);
    

    in.close();

}
#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check;

void help()
{
	cout << "k_vs_time." << endl;

	cout << "Usage: k_vs_time [options] <index> <patterns>" << endl;
	cout << "   <index>      index file (with extension .bri)" << endl;
	cout << "   <pattern>    test pattern file" << endl;

	exit(0);
}


template<class T>
void locate_all(ifstream& in, string patterns)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    auto t1 = high_resolution_clock::now();

    T idx;

    cout << "Loading br-index ..." << endl;

    idx.load(in);

    cout << "done." << endl;

    auto t2 = high_resolution_clock::now();


    ifstream ifs(patterns);

    //read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
    string header;
    getline(ifs,header);

    ulint n = get_number_of_patterns(header);
    ulint m = get_patterns_length(header);

    cout << "patterns: " << n << " pattern length: " << m << endl;

    ulint last_perc = 0;

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    auto t5 = high_resolution_clock::now();
    ulint count_time = 0;
    ulint locate_time = 0;
    ulint tot_time = 0;


    vector<ulint> times(m+1,0);

    for (ulint i = 0; i < n; ++i)
    {
        ulint perc = 100 * i / n;
        if (perc > last_perc)
        {
            cout << perc << "% done ..." << endl;
            last_perc = perc;
        }

        string p = string();

        for (ulint j = 0; j < m; ++j)
        {
            char c;
            ifs.get(c);
            p += c;
        }

        br_sample sample(idx.get_initial_sample());
        for (ulint j = 0; j < m; ++j)
        {
            sample = idx.right_extension(p[j],sample);
        }

        if (sample.is_invalid())
        {
            cout << "Invalid pattern. Go next." << endl; continue;
        }

        for (ulint j = m; j >= 2; --j)
        {
            t3 = high_resolution_clock::now();
            sample = idx.left_contraction(sample);
            t4 = high_resolution_clock::now();
            times[j] += duration_cast<microseconds>(t4-t3).count();
        }
        if (sample.len != 1) cout << "sample.len is not 1. error." << endl;

    }



    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;
    cout << "#Patterns      n = " << n << endl;
    cout << "Pattern length m = " << m << endl;

    cout << "Pattern contraction times from k=2 to k=m:" << endl;
    cout << "[";
    for (ulint k = 2; k < m; ++k)
    {
        cout << times[k] << ", ";
    }
    cout << times[m] << "]" << endl;
    
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    string idx_file(argv[1]);
    string patt_file(argv[2]);



    ifstream in(idx_file);

    cout << "Loading br-index" << endl;

    locate_all<br_index<> >(in, patt_file);
    
    

    in.close();

}
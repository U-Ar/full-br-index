#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index_fixed.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check;

void help()
{
	cout << "fixedtest." << endl;

	cout << "Usage: bri-mem [options] <index> <patterns>" << endl;
    cout << "   <textfile>   original text file" << endl;
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

    idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "testing de Bruijn version ..." << endl;
    

    ifstream ifs(patterns);

    stringstream ssp;
    ssp << ifs.rdbuf();
    string p = ssp.str();


    ulint n = 1;
    ulint m = p.size();

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    auto t5 = high_resolution_clock::now();
    ulint count_time = 0;
    ulint locate_time = 0;
    ulint tot_time = 0;



    t3 = high_resolution_clock::now();

    ulint errors = 0;
    ulint length = 64;
    ulint iterations = 1000;

    br_sample fixed_sample(idx.get_initial_sample());
    for (ulint i = 0; i < iterations; ++i)
    {
        br_sample sample(idx.get_initial_sample());
        for (ulint j = 0; j < length; ++j)
        {
            fixed_sample = idx.outgoing(p[i*length+j],fixed_sample);
            sample = idx.right_extension(p[i*length+j],sample);
        }
        cout << i << "-th iteration:  fixed occ: " << fixed_sample.size() << " normal occ: " << sample.size() << endl;
        if (fixed_sample.size() != sample.size()) 
        {
            errors++; continue;
        }

        auto fixed_occs = idx.locate_sample(fixed_sample);
        auto occs = idx.locate_sample(sample);
        std::sort(fixed_occs.begin(),fixed_occs.end());
        std::sort(occs.begin(),occs.end());
        for (ulint k = 0; k < occs.size(); ++k)
        {
            if (occs[k] != fixed_occs[k]) errors++;
        }

        cout << "pattern: ";
        for (ulint l = 0; l < length; ++l) cout << p[i*length+l];
        cout << endl;

        string label = idx.label(fixed_sample);

        for (ulint k = 0; k < occs.size(); ++k)
        {
            ulint occ = occs[k];
            for (ulint l = 0; l < length; ++l)
            {
                if (check[occ+l] != p[i*length+l]) errors++;
                if (check[occ+l] != label[l]) errors++;
            }
        }

    }

    fixed_sample = idx.get_initial_sample();
    for (ulint i = iterations; i-- > 0;)
    {
        br_sample sample(idx.get_initial_sample());
        for (ulint j = length; j-- > 0;)
        {
            fixed_sample = idx.incoming(p[i*length+j],fixed_sample);
            sample = idx.left_extension(p[i*length+j],sample);
        }
        cout << iterations-i << "-th iteration:  fixed occ: " << fixed_sample.size() << " normal occ: " << sample.size() << endl;
        if (fixed_sample.size() != sample.size())
        {
            errors++; continue;
        }


        auto fixed_occs = idx.locate_sample(fixed_sample);
        auto occs = idx.locate_sample(sample);
        std::sort(fixed_occs.begin(),fixed_occs.end());
        std::sort(occs.begin(),occs.end());
        for (ulint k = 0; k < occs.size(); ++k)
        {
            if (occs[k] != fixed_occs[k]) errors++;
        }

        cout << "pattern: ";
        for (ulint l = 0; l < length; ++l) cout << p[i*length+l];
        cout << endl;

        string label = idx.label(fixed_sample);

        for (ulint k = 0; k < occs.size(); ++k)
        {
            ulint occ = occs[k];
            for (ulint l = 0; l < length; ++l)
            {
                if (check[occ+l] != p[i*length+l]) errors++;
                if (check[occ+l] != label[l]) errors++;
            }
        }
    }


    cout << "Number of errors: " << errors << endl;



    t4 = high_resolution_clock::now();


    t5 = high_resolution_clock::now();

    count_time += duration_cast<microseconds>(t4-t3).count();
    locate_time += duration_cast<microseconds>(t5-t4).count();
    tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();


    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    cout << "Pattern length     m =  " << m << endl;
    
    cout << "LF-mapping time: " << count_time << " microseconds" << endl;
    cout << "Phi        time: " << locate_time << " microseconds" << endl;
    cout << "Total time : " << tot_time << " microseconds" << endl;
    
}



int main(int argc, char** argv)
{
    if (argc < 4) help();

    string check_file(argv[1]);
    string idx_file(argv[2]);
    string patt_file(argv[3]);


    ifstream in_check(check_file);

    stringstream ss_check;
    ss_check << in_check.rdbuf();
    check = ss_check.str();

    ifstream in(idx_file);

    cout << "Loading br-index-fixed" << endl;

    locate_all<br_index_fixed<> >(in, patt_file);
    
    

    in.close();

}
#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "br_index_limited.hpp"
#include "br_index_st.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check = string();
bool limited = false;
bool once = false;

void help()
{
	cout << "bri-mem: locate maximal exact matches of the input patterns. " << endl;
    cout << "         by default br_index_st is used." << endl;

	cout << "Usage: bri-mem [options] <index> <patterns>" << endl;
    cout << "   -limited     use the limited version of br-index." << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file" << endl;
    cout << "   -once        use a text file as one pattern instead of pattern file" << endl;
	cout << "   <index>      index file (with extension .bri)" << endl;
	cout << "   <patterns>   file in pizza&chili format or simple text file(if -once) containing the patterns." << endl;

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
    else if (s.compare("-once") == 0)
    {
        once = true;
    }
    else if (s.compare("-limited") == 0)
    {
        limited = true;
    }
    else
    {

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

template<class T>
void locate_all(ifstream& in, string patterns)
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

    auto t1 = high_resolution_clock::now();

    T idx;

    idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "searching maximal exact matches ..." << endl;
    
    if (!once)
    {
        ifstream ifs(patterns);

        //read header of the pizza&chilli input file
        //header example:
        //# number=7 length=10 file=genome.fasta forbidden=\n\t
        string header;
        getline(ifs,header);

        ulint n = get_number_of_patterns(header);
        ulint m = get_patterns_length(header);

        ulint last_perc = 0;

        ulint occ_tot = 0;

        auto t3 = high_resolution_clock::now();
        auto t4 = high_resolution_clock::now();
        auto t5 = high_resolution_clock::now();
        ulint count_time = 0;
        ulint locate_time = 0;
        ulint tot_time = 0;

        // extract patterns from file and search them in the index
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

            t3 = high_resolution_clock::now();

            idx.maximal_exact_match(p);

            t4 = high_resolution_clock::now();


            t5 = high_resolution_clock::now();

            count_time += duration_cast<microseconds>(t4-t3).count();
            locate_time += duration_cast<microseconds>(t5-t4).count();
            tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();

        }
        
        
        ifs.close();

        ulint load = duration_cast<milliseconds>(t2-t1).count();
        cout << "Load time  : " << load << " milliseconds" << endl;

        cout << "Number of patterns n = " << n << endl;
        cout << "Pattern length     m = " << m << endl;
        
        cout << "LF-mapping time: " << count_time << " microseconds" << endl;
        cout << "Phi        time: " << locate_time << " microseconds" << endl;
        cout << "Total time : " << tot_time << " microseconds" << endl;
        cout << "Search time: " << (double)tot_time/n << " microseconds/pattern (total: " << n << " patterns)" << endl;
    }
    else if (once)
    {

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

        /*std::pair<
            ulint,
            std::vector<std::pair<ulint,br_sample>>
        > res = idx.maximal_exact_match(p);*/
        idx.maximal_exact_match(p);

        t4 = high_resolution_clock::now();

        /*std::vector<ulint> occs;
        for (auto p: res.second) {
            auto tmp = idx.locate_sample(p.second);
            occs.insert(occs.end(),tmp.begin(),tmp.end());
        }*/

        t5 = high_resolution_clock::now();

        count_time += duration_cast<microseconds>(t4-t3).count();
        locate_time += duration_cast<microseconds>(t5-t4).count();
        tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();


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

        cout << "Pattern length     m =  " << m << endl;
        
        cout << "LF-mapping time: " << count_time << " microseconds" << endl;
        cout << "Phi        time: " << locate_time << " microseconds" << endl;
        cout << "Total time : " << tot_time << " microseconds" << endl;
    }
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

    cout << "Loading br-index" << endl;

    if (limited)
    {
        locate_all<br_index_limited<> >(in, patt_file);
    }
    else 
    {
        locate_all<br_index_st<> >(in, patt_file);
    }
    

    in.close();

}
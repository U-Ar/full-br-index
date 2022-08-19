#include <iostream>
#include <chrono>
#include <cstdlib>
#include <random>
#include <string>

#include "br_index.hpp"
#include "br_index_limited.hpp"
#include "br_index_st.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check = string();
string text;

void help()
{
	cout << "x_treeop: test suffix tree operations." << endl;

	cout << "Usage: x_treeop -c <text> <index>" << endl;
    cout << "   -c <text>    original text" << endl;
    cout << "                used to randomly choose ST nodes" << endl;
    cout << "   <index>      index file of r-ST" << endl;

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
    else
    {

		cout << "Error: unknown option " << s << endl;
		help();

	}
}



void tree_op_1(br_index_st<>& idx, std::mt19937& engine)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;
    using std::chrono::nanoseconds;

    std::uniform_int_distribution<> dist(0,idx.text_size()-1);

    long long int traversed = 0;
    long long int sdepth_time = 0;
    long long int child_time = 0;
    long long int fchild_time = 0;
    long long int nsibling_time = 0;
    long long int parent_time = 0;

    while (traversed < 10000)
    {
        // choose a random leaf
        int pos = dist(engine);
        auto sample = idx.get_initial_sample();
        for (ulint k = pos; k < idx.text_size(); ++k)
        {
            sample = idx.right_extension(text[k],sample);
            if (sample.is_leaf()) break;
        }
        // std::cout << "pos:  " << pos << std::endl;
        // std::cout << "range:" << sample.range.first << "," << sample.range.second << " len:" << sample.len << std::endl;
        // std::cout << "j:" << sample.j << "  d:" << sample.d << std::endl;

        assert(sample.is_leaf());
        assert(!sample.is_invalid());

        while (sample.len != 0)
        {
            // sdepth, child, fchild, nsibling, parent
            auto t1 = high_resolution_clock::now();
            auto depth = idx.sdepth(sample);
            auto t2 = high_resolution_clock::now();
            assert(depth <= idx.bwt_size());

            //auto t3 = high_resolution_clock::now();
            //auto csample = idx.child(sample,text[dist(engine)]);
            //auto t4 = high_resolution_clock::now();
            //assert(csample.range.second <= idx.bwt_size());

            auto t5 = high_resolution_clock::now();
            auto fcsample = idx.fchild(sample);
            auto t6 = high_resolution_clock::now();
            assert(fcsample.range.second <= idx.bwt_size());

            auto t7 = high_resolution_clock::now();
            auto nssample = idx.nsibling(sample);
            auto t8 = high_resolution_clock::now();
            assert(nssample.range.second <= idx.bwt_size());

            auto t9 = high_resolution_clock::now();
            sample = idx.parent(sample);
            auto t10 = high_resolution_clock::now();
            assert(!sample.is_invalid());

            traversed++;
            sdepth_time += duration_cast<nanoseconds>(t2-t1).count();
            //child_time += duration_cast<nanoseconds>(t4-t3).count();
            fchild_time += duration_cast<nanoseconds>(t6-t5).count();
            nsibling_time += duration_cast<nanoseconds>(t8-t7).count();
            parent_time += duration_cast<nanoseconds>(t10-t9).count();
            if (traversed > 10000 || sample.len == 0) break;
        }
    }
    cout << "(" << traversed << " traversed times)" << endl;
    cout << "sdepth  : " << (double)sdepth_time / traversed / 1000 << " microseconds" << endl;
    //cout << "child   : " << (double)child_time / traversed / 1000 << " microseconds" << endl;
    cout << "fchild  : " << (double)fchild_time / traversed / 1000 << " microseconds" << endl;
    cout << "nsibling: " << (double)nsibling_time / traversed / 1000 << " microseconds" << endl;
    cout << "parent  : " << (double)parent_time / traversed / 1000 << " microseconds" << endl;
}

void tree_op_2(br_index_st<>& idx, std::mt19937& engine)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;
    using std::chrono::nanoseconds;

    std::uniform_int_distribution<> dist(0,idx.text_size()-1);

    int traversed = 0;
    long long int letter_time = 0;

    while (traversed < 10000)
    {
        // choose a random leaf
        int pos = dist(engine);
        auto sample = idx.get_initial_sample();
        for (ulint k = pos; k < idx.text_size(); ++k)
        {
            sample = idx.right_extension(text[k],sample);
            if (sample.is_leaf()) break;
        }

        assert(!sample.is_invalid());

        while (sample.len >= 5)
        {
            // letter
            auto t1 = high_resolution_clock::now();
            auto character = idx.letter(sample,4);
            auto t2 = high_resolution_clock::now();
            assert(character == text[sample.j-sample.d+4]);

            sample = idx.parent(sample);
            traversed++;
            letter_time += duration_cast<nanoseconds>(t2-t1).count();
            if (traversed > 10000) break;
        }
    }
    cout << "(" << traversed << " traversed times)" << endl;
    cout << "letter: " << (double)letter_time / traversed / 1000 << " microseconds" << endl;
}

void tree_op_3(br_index_st<>& idx, std::mt19937& engine)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;
    using std::chrono::nanoseconds;

    std::uniform_int_distribution<> dist(0,idx.text_size()-1);

    int traversed = 0;
    long long int total_time = 0;

    while (traversed < 10000)
    {
        // choose a random leaf
        int pos = dist(engine);
        auto sample = idx.get_initial_sample();
        for (ulint k = pos; k < idx.text_size(); ++k)
        {
            sample = idx.right_extension(text[k],sample);
            if (sample.is_leaf()) break;
        }
        sample = idx.parent(sample);

        assert(!sample.is_invalid());

        while (sample.len != 0)
        {
            auto t1 = high_resolution_clock::now();
            sample = idx.slink(sample);
            auto t2 = high_resolution_clock::now();
            assert(!sample.is_invalid());

            traversed++;
            total_time += duration_cast<nanoseconds>(t2-t1).count();
            if (traversed > 10000) break;
        }
    }
    cout << "(" << traversed << " traversed times)" << endl;
    cout << "suffix-link: " << (double)total_time / traversed / 1000 << " microseconds" << endl;
}

void tree_op_4(br_index_st<>& idx, std::mt19937& engine)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;
    using std::chrono::nanoseconds;

    std::uniform_int_distribution<> dist(0,idx.text_size()-1);

    long long int total_time = 0;

    for (int traversed = 0; traversed < 10000; ++traversed)
    {
        int pos1 = dist(engine);
        auto sample1 = idx.get_initial_sample();
        for (ulint k = pos1; k < idx.text_size(); ++k)
        {
            sample1 = idx.right_extension(text[k],sample1);
            if (sample1.is_leaf()) break;
        }
        assert(sample1.is_leaf());
        assert(!sample1.is_invalid());


        int pos2 = dist(engine);
        auto sample2 = idx.get_initial_sample();
        for (ulint k = pos2; k < idx.text_size(); ++k)
        {
            sample2 = idx.right_extension(text[k],sample2);
            if (sample2.is_leaf()) break;
        }
        assert(sample2.is_leaf());
        assert(!sample2.is_invalid());


        auto t1 = high_resolution_clock::now();
        auto sample_lca = idx.lca(sample1,sample2);
        auto t2 = high_resolution_clock::now();
        assert(!sample_lca.is_invalid());
        sample_lca.range.first = sample_lca.range.second;

        total_time += duration_cast<nanoseconds>(t2-t1).count();
    }

    // per op /10000 & nano->micro /1000
    cout << "LCA: " << (double)total_time / 10000000 << " microseconds" << endl;
}

void tree_op_5(br_index_st<>& idx, std::mt19937& engine)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;
    using std::chrono::nanoseconds;

    std::uniform_int_distribution<> dist(0,idx.text_size()-1);

    int traversed = 0;
    long long int child_time = 0;
    long long int bchild_time = 0;

    while (traversed < 10000)
    {
        // choose a random leaf
        int pos = dist(engine);
        auto sample = idx.get_initial_sample();
        for (ulint k = pos; k < idx.text_size(); ++k)
        {
            sample = idx.right_extension(text[k],sample);
            if (sample.is_leaf()) break;
        }
        sample = idx.parent(sample);

        assert(!sample.is_invalid());

        while (sample.len != 0)
        {
            // sdepth, child, fchild, nsibling, parent
            auto children = idx.child_chars(sample);
            int degree = children.size();
            if (degree >= 3)
            {
                std::uniform_int_distribution<> char_dist(0,degree-1);
                unsigned char ch = children[char_dist(engine)];

                auto t1 = high_resolution_clock::now();
                auto child = idx.child(sample,ch);
                auto t2 = high_resolution_clock::now();
                assert(!child.is_invalid());

                traversed++;
                child_time += duration_cast<nanoseconds>(t2-t1).count();
            }
            

            sample = idx.parent(sample);
            
            if (traversed > 10000) break;
        }
    }
    cout << "(" << traversed << " traversed times)" << endl;
    cout << "child     : " << (double)child_time / traversed / 1000 << " microseconds" << endl;
}


void tree_op(br_index_st<>& idx)
{
    cout << "tree operations on br_index_st." << endl;
    cout << "1. from a random leaf to the root, 10000 traversals" << endl;
    cout << "    1.1 sdepth" << endl;
    cout << "    1.2 fchild" << endl;
    cout << "    1.3 nsibling" << endl;
    cout << "    1.4 parent" << endl;
    cout << "2. same as 1 but depth>=5 is only allowed" << endl;
    cout << "    2.1 letter" << endl;
    cout << "3. from a parent of a random leaf, 10000 suffix links" << endl;
    cout << "    3.1 slink" << endl;
    cout << "    3.2 tdepth" << endl;
    cout << "4. 10000 random leaf pairs" << endl;
    cout << "    4.1 lca" << endl << endl;
    cout << "5. only deg >= 3 && selecting random child" << endl;
    cout << "    5.1 child" << endl << endl;

    std::seed_seq seed{2022};
    std::mt19937 engine(seed);
    // std::uniform_int_distribution<> dist(0,9);
    // dist(engine)

    // 1. from a random leaf to the root
    tree_op_1(idx,engine);
    // 2. same as 1 but only sdepth>=5 is allowed
    tree_op_2(idx,engine);
    // 3. from a parent of a random leaf, 10000 suffix links
    tree_op_3(idx,engine);
    // 4. 10000 random leaf pairs
    tree_op_4(idx,engine);
    // 5. child
    tree_op_5(idx,engine);

    
}


int main(int argc, char** argv)
{
    if (argc < 2) help();

    int ptr = 1;

    while (ptr < argc-1) parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);

    cout << "Loading br-index ... " << flush;
    ifstream in(idx_file);
    br_index_st<> idx;
    idx.load(in);
    cout << "done." << endl;

    bool c = false;
    if (check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    tree_op(idx);
    
}
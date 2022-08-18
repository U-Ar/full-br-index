#include<iostream>
#include <chrono>

#include "br_index_fixed.hpp"

using namespace std;
using namespace bri;

long length = 0;

void help(){
	cout << "window_slide: " << endl;
	cout << "Usage:       window_slide <text> <index>" << endl;
	cout << "   <text>    original text" << endl;
	cout << "   <index>   used index" << endl;
    cout << "   <length>  max window length" << endl;
	exit(0);
}


/*

Experiment
If we execute this on length 1 index, naive method is applied
If we execute this on length 16 index, shortcut is applied until l=16 window

- naive: from_empty
- bi:    right_contraction -> left_extension
- uni:   right_contraction_uni -> left_extension_uni

*/


int main(int argc, char** argv){
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

	if(argc < 4)
		help();

    
    ifstream ifs(argv[1]);
    stringstream ss;
    ss << ifs.rdbuf();
    string text = ss.str();

	char* e;
	length = strtol(argv[3],&e,10);
	if(*e != '\0' || length < 0){
		cout << "Error: invalid fixed length." << endl;
		return 0;
	}


	br_index_fixed<> idx;
	cout << "Loading br-index-fixed from" << argv[2] << endl;
	idx.load_from_file(argv[2]);
    cout << "index border length: " << idx.border_length() << endl;

    auto ts_naive = high_resolution_clock::now();
    auto te_naive = high_resolution_clock::now();
    auto ts_bi = high_resolution_clock::now();
    auto te_bi = high_resolution_clock::now();
    auto ts_uni = high_resolution_clock::now();
    auto te_uni = high_resolution_clock::now();

    cout << endl << "Experiments on l=2 to " << length << endl;
	for (ulint l = 2; l <= length; ++l)
    {
        cout << "l=" << l << endl;

        ulint time_naive = 0;
        ulint time_bi = 0;
        ulint time_uni = 0;

        // prepare samples
        // auto sample_naive = idx.get_initial_sample();
        auto sample_bi = idx.get_initial_sample();
        auto sample_uni = idx.get_initial_sample();
        for (ulint i = text.size()-1; i >= text.size()-l; --i)
        {
            sample_bi = idx.left_extension(text[i],sample_bi);
            sample_uni = idx.left_extension_uni(text[i],sample_uni);
        }
        assert(sample_bi.len == l && sample_uni.len == l && sample_bi.size()>0 && sample_uni.size()>0);

        for (ulint i = text.size()-l; i-- > 0;)
        {
            ts_naive = high_resolution_clock::now();
            auto sample_naive = idx.get_initial_sample();
            for (ulint d = 0; d < l; ++d) {
                sample_naive = idx.left_extension_uni(text[i+l-1-d],sample_naive);
            }
            te_naive = high_resolution_clock::now();
            time_naive += duration_cast<milliseconds>(te_naive-ts_naive).count();


            ts_bi = high_resolution_clock::now();
            sample_bi = idx.right_contraction(sample_bi);
            sample_bi = idx.left_extension(text[i],sample_bi);
            te_bi = high_resolution_clock::now();
            time_bi += duration_cast<milliseconds>(te_bi-ts_bi).count();


            ts_uni = high_resolution_clock::now();
            sample_uni = idx.right_contraction_uni(sample_uni);
            sample_uni = idx.left_extension_uni(text[i],sample_uni);
            te_uni = high_resolution_clock::now();
            time_uni += duration_cast<milliseconds>(te_uni-ts_uni).count();


            assert(sample_naive.range.first == sample_bi.range.first);
            assert(sample_naive.range.first == sample_uni.range.first);
            assert(sample_naive.range.second == sample_bi.range.second);
            assert(sample_naive.range.second == sample_uni.range.second);
        }

        std::cout << "<naive> " << time_naive << std::endl;
        std::cout << "<bi>    " << time_bi    << std::endl;
        std::cout << "<uni>   " << time_uni   << std::endl;
    }

}
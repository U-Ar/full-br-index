#include "br_index.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>

using namespace bri;
using namespace std;

const ulint LENGTH = 50;

struct KR_window {
    int wsize;
    //int *window;
    vector<int> window;
    int asize;
    // const ulint prime = 1999999973;
    const ulint prime = 27162335252586509;
    ulint hash;
    ulint tot_char;
    ulint asize_pot;   // asize^(wsize-1) mod prime
    
    KR_window(int w): wsize(w) {
        asize = 256;
        asize_pot = 1;
        for(int i=1;i<wsize;i++) 
        asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm  
        // alloc and clear window
        // window = new int[wsize];
        window = vector<int>(w);
        reset();     
    }

    void resize(int w) {
        wsize = w;
        asize = 256;
        asize_pot = 1;
        for(int i=1;i<wsize;i++) 
        asize_pot = (asize_pot*asize)% prime;
        // delete[] window;
        // window = new int[wsize];
        window.resize(w);
        reset();
    }
    
    // init window, hash, and tot_char 
    void reset() {
        for(int i=0;i<wsize;i++) window[i]=0;
        // init hash value and related values
        hash=tot_char=0;    
    }
    
    ulint addchar(int c) {
        int k = tot_char++ % wsize;
        // complex expression to avoid negative numbers 
        hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution  
        hash = (asize*hash + c) % prime;      //  add char i 
        window[k]=c;
        // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
        return hash; 
    }
    // debug only 
    string get_window() {
        string w = "";
        int k = (tot_char-1) % wsize;
        for(int i=k+1;i<k+1+wsize;i++)
        w.append(1,window[i%wsize]);
        return w;
    }

};

ulint kr_hash(string const& s) {
    ulint hash = 0;
    //const ulint prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const ulint prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k = 0; k < s.size(); k++) {
        int c = (unsigned char) s[k];
        assert(c>=0 && c< 256);
        hash = (256*hash + c) % prime;    //  add char k
    } 
    return hash; 
}

// get frequencies of substrings of length <= max_len 
// keys are Karp-Rabin hash
map<ulint,ulint> kr_freqs(string const& text, ulint max_len)
{
    vector<KR_window> windows(max_len, KR_window(1));
    for (ulint l = 1; l <= max_len; ++l) windows[l-1].resize(l);

    map<ulint,ulint> hash_table;

    for (size_t i = 0; i < text.size(); ++i)
    {
        for (ulint l = 1; l <= max_len; ++l) 
        {
            ulint hash = windows[l-1].addchar((uchar)text[i]); 
            if (i >= l-1) {
                // increment hash table's freq
                if (hash_table.count(hash) == 0) {
                    hash_table[hash] = 1;
                } else {
                    hash_table[hash]++;
                }
            }
        }
    }

    return hash_table;
}

bool verify(string const& text, br_index& idx, br_sample const& sample,
    string const& pattern, map<ulint, ulint>& freqs,
    ulint i, ulint l, bool rightward)
{
    bool ok = true;
    string opname = rightward ? "rightward MEMs" : "leftward MEMs";

    ulint cnt = idx.count_sample(sample);

    if (l <= LENGTH) {
        ulint hash = kr_hash(pattern.substr(i,l));

        // check frequency of substring
        if (freqs.find(hash)==freqs.end() || freqs.at(hash)!=cnt)
        {
            cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
            cerr << "  Numbers of substring \"" << pattern.substr(i,l) << "\" differ." << endl;
            cerr << "  KR hash  : " << (freqs.find(hash)==freqs.end() ? 0 : freqs.at(hash)) << endl;
            cerr << "  br-index : " << cnt;
            return false;
        }
    }

    // check non-existence of left-extension
    if (i > 0) {
        br_sample ex_sample(idx.left_extension(sample,(uchar)pattern[i-1]));
        if (ex_sample.is_valid()) {
            cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
            cerr << "  current pattern: " << pattern.substr(i,l) << endl;
            cerr << "    is not maximal." << endl;
            cerr << "  occurring character on the left: " << pattern[i-1] << endl;
            cerr << "  number of occurrence           : " << ex_sample.size() << endl;
            return false;
        }
    }

    // check non-existence of right-extension
    if (i+l < pattern.size()) {
        br_sample ex_sample(idx.right_extension(sample,(uchar)pattern[i+l]));
        if (ex_sample.is_valid()) {
            cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
            cerr << "  current pattern: " << pattern.substr(i,l) << endl;
            cerr << "    is not maximal." << endl;
            cerr << "  occurring character on the right: " << pattern[i+l] << endl;
            cerr << "  number of occurrence            : " << ex_sample.size() << endl;
            return false;
        }
    }

    // check locate correctness
    vector<ulint> positions(idx.locate_sample(sample));
    if (positions.size() != cnt) {
        cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
        cerr << "  Size of location vector differ from count val." << endl;
        return false;
    }

    for (ulint pos : positions) {
        // check text position overflow
        if (pos+l>text.size()) {
            cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
            cerr << "  Pattern starting from located pos overflows text size." << endl;
            cerr << "  pos: " << pos << " pattern: " << pattern.substr(i,l) << " textsize: " << text.size() << endl;
            ok = false; continue;
        }

        // check the match between the pattern and located substrings
        string reference = text.substr(pos,l);
        if (pattern.substr(i,l) != reference) {
            cerr << "Error at " << opname << ", position " << i << ", length " << l << "." << endl;
            cerr << "  Substring starting from located pos differ from the pattern." << endl;
            cerr << "  Pos : " << pos << endl;
            cerr << "  Substring at pos        : " << reference << endl;
            cerr << "  Pattern we are searching: " << pattern.substr(i,l) << endl;
            ok = false;
        }
    }

    return ok;
}

int main(int argc, char** argv)
{
    using clock = chrono::high_resolution_clock;

    auto t1 = clock::now();

    if (argc != 5) { cerr << "Invalid command line argument. Exitting ... " << endl; exit(1); }

    const string EXTIDX = ".bri";
    const string EXTIN = ".in";
    const string EXTPFP = ".pfp";

    string input_file = argv[1];
    string idx_file = argv[2];
    string patt_file = argv[3];
    long bl = stoi(argv[4]);

    cout << "==== Test2 on " << input_file;
    if (idx_file.find(EXTIN+EXTIDX) != string::npos) cout << " with in-memory construction, ";
    else cout << " with PFP construction, ";
    cout << "bl=" << bl << ", pattern file: " << patt_file << endl;

    // read input file
    ifstream fin(input_file);
    if (!fin.is_open()) {
        cerr << "Text file " << input_file << " is not found. Exitting ... " << endl; exit(1);
    }

    string text((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());
    fin.close();


    // read pattern
    ifstream fp(patt_file);
    if (!fp.is_open()) {
        cerr << "Pattern file " << patt_file << " is not found. Exitting ... " << endl; exit(1);
    }
    string pattern;
    while (!fp.eof()) pattern.push_back(fp.get());
    fp.close();


    // compute KR hash
    cout << "Count frequencies of substrings of length <= " << LENGTH << endl;
    map<ulint,ulint> freqs = kr_freqs(text, LENGTH);


    // read index file
    ifstream fidx(idx_file);
    if (!fidx.is_open()) {
        cerr << "Index file " << idx_file << " is not found. Exitting ... " << endl; exit(1);
    }
    cout << "Loading " << idx_file << " with bl = " << bl << " ... " << flush;
    br_index idx;
    idx.load(fidx,bl);
    fidx.close();
    cout << "done." << endl;


    if (idx.text_size() != text.size()) {
        cerr << "Text size and index's text size differ." << endl;
        cerr << "Text size:         " << text.size() << endl;
        cerr << "Index's text size: " << idx.text_size() << "  Exitting ..." << endl;
        exit(1);
    }

    // count index size
    ulint idx_size = idx.get_space();
    cout << "Index size: " << idx_size << " bytes" << endl;

    ulint n = text.size();
    ulint m = pattern.size();

    bool all_ok = true;

    // compute MEMs from left to right
    {
        ulint j = 0, l = 0;
        bool extended = false, ok = true;

        br_sample sample(idx.get_initial_sample());
        
        for (ulint i = 0; i < m; ++i) {
            while (j < m)
            {
                br_sample new_sample = idx.right_extension(sample, (uchar)pattern[j]);
                if (new_sample.is_invalid()) break;
                extended = true;
                sample = new_sample;
                l++;
                j++;
            }
            if (extended) // maximal substring
            {
                ok = verify(text, idx, sample, pattern, freqs, i, l, true);
            }
            if (!ok) break;
            
            if (i == j)
            {
                sample = idx.get_initial_sample();
                j++;
                l = 0;
            }
            else 
            {
                sample = idx.left_contraction(sample);
                l--;
            }
            extended = false;
        }
        all_ok = ok && all_ok;
    }

    // compute MEMs from right to left
    {
        long j = m-1;
        ulint l = 0;
        bool extended = false, ok = true;

        br_sample sample(idx.get_initial_sample());
        
        for (ulint i = m; i-- > 0;) {
            while (j >= 0)
            {
                br_sample new_sample = idx.left_extension(sample, (uchar)pattern[j]);
                if (new_sample.is_invalid()) break;
                extended = true;
                sample = new_sample;
                l++;
                j--;
            }
            if (extended) // maximal substring
            {
                ok = verify(text, idx, sample, pattern, freqs, i-l+1, l, false);
            }
            if (!ok) break;
            
            if (i == (ulint)j)
            {
                sample = idx.get_initial_sample();
                j--;
                l = 0;
            }
            else 
            {
                sample = idx.right_contraction(sample);
                l--;
            }
            extended = false;
        }
        all_ok = ok && all_ok;
    }

    auto t2 = clock::now();

    cout << "==== Test2 on " << input_file << " finished" << endl;
    cout << "==== (Elapsed time: " 
         << fixed << setprecision(4) 
         << (double)chrono::duration_cast<chrono::milliseconds>(t2-t1).count()/1000
         << ")" << endl << endl;

    if (!all_ok) exit(1);
}
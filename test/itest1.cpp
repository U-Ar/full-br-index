#include "br_index.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

using namespace bri;
using namespace std;

struct KR_window {
  int wsize;
  int *window;
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
    window = new int[wsize];
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
  
  ~KR_window() {
    delete[] window;
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
    vector<KR_window> windows(max_len);
    for (ulint l = 1; l <= max_len; ++l) windows[l-1] = KR_window(l);

    map<ulint,ulint> hash_table;

    for (size_t i = 0; i < text.size(); ++i)
    {
        for (ulint l = 1; l <= max_len; ++l) 
        {
            ulint hash = windows[l-1].addchar(text[i]); 
            if (i >= l-1) {
                // increment hash table's freq
                if (hash_table.find(hash)==hash_table.end()) {
                    hash_table[hash] = 1;
                } else {
                    hash_table[hash]++;
                }
            }
        }
    }

    return hash_table;
}

int main(int argc, char** argv)
{
    if (arc != 4) { cerr << "Invalid command line argument. Exitting ... " << endl; exit(1); }

    const ulint LENGTH = 50;

    const string EXTIDX = ".bri";
    const string EXTIN = ".in";
    const string EXTPFP = ".pfp";

    string input_file = argv[1];
    string idx_file = argv[2];
    long bl = stoi(argv[3]);

    cout << "==== Test1 on " << input_file;
    if (idx_file.find(EXTIN+EXTIDX) != string::npos) cout << " with in-memory construction, ";
    else cout << " with PFP construction, ";
    cout << "bl=" << bl << endl;


    // read input file
    ifstream fin(input_file);
    if (!fin.is_open()) {
        cerr << "Text file " << input_file << " is not found. Exitting ... " << endl; exit(1);
    }
    string text;
    while (!fin.eof()) text.push_back(fin.get());
    fin.close();


    // compute KR hash
    cout << "Count frequencies of substrings of length <= " << LENGTH << endl;
    map<ulint,ulint> freqs = kr_freqs(text, LENGTH);
    

    // read index file
    ifstream fidx(idx_file);
    if (!fidx.is_open()) {
        cerr << "Index file " << idx_file << " is not found. Exitting ... " << endl; exit(1);
    }
    cout << "Loading " << idx_file << " with bl = " << bl << endl;
    br_index idx;
    idx.load(fidx,bl);

    if (idx.text_size() != text.size()) {
        cerr << "Text size and index's text size differ." << endl;
        cerr << "Text size:         " << text.size() << endl;
        cerr << "Index's text size: " << idx.text_size() << "  Exitting ..." << endl;
        exit(1);
    }

    // count index size
    ulint idx_size = idx.get_space();
    cout << "Index size: " << idx_size << " bytes" << endl;


    ulint fragments = text.size() / LENGTH;
    cout << "Divide texts into " << fragments << " substrings" << endl;
    cout << "Substring length: " << LENGTH << endl;

    // test for each fragment
    for (ulint k = 0; k < fragments; ++k) {
        string pattern = text.substr(k*LENGTH,LENGTH);

    }












}
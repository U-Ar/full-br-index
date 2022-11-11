#include "br_index.hpp"

#include <cstdlib>
#include <iostream>

using namespace bri;
using namespace std;

const ulint LENGTH = 50;

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

int main(int argc, char** argv)
{
    if (arc != 4) { cout << "Invalid command line argument." << endl; exit(1); }

    string input_file = argv[1];
    string idx_file = argv[2];
    long bl = stoi(argv[3]);


}
/*
 * fully functional bidirectional r-index
 * author: Yuma Arakawa
 */

#ifndef INCLUDED_BR_INDEX_FULL_HPP
#define INCLUDED_BR_INDEX_FULL_HPP

#include "definitions.hpp"
#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include "permuted_lcp.hpp"

namespace bri {

class br_index {
    
public:
    using triple = std::tuple<range_t, ulint, ulint>;
    using sparse_bitvector_t = sparse_sd_vector;
    using rle_string_t = rle_string<>;

    friend class br_index_builder;


    br_index() {}
    br_index(std::string const& input, ulint length=8, bool sais = true);

    range_t full_range();
    ulint LF(ulint i);
    range_t LF(range_t rn, uchar c); // c: remapped char
    ulint LFR(ulint i);
    range_t LFR(range_t rn, uchar c); // c: remapped char
    ulint FL(ulint i);
    ulint FLR(ulint i);
    ulint Phi(ulint i);
    ulint PhiR(ulint i);
    ulint PhiI(ulint i);
    ulint PhiIR(ulint i);

    uchar F_at(ulint i);
    range_t get_char_range(uchar c); // c: remapped char


    // search initialization
    br_sample get_initial_sample();
    inline br_sample get_invalid_sample();
    inline bool is_invalid(br_sample const& sample) const { return sample.is_invalid(); }


    // br-index functionality
    br_sample left_extension(br_sample const& prev_sample, uchar c);
    br_sample right_extension(br_sample const& prev_sample, uchar c);
    br_sample left_contraction(br_sample const& prev_sample);
    br_sample right_contraction(br_sample const& prev_sample);

    br_sample search(std::string const& pattern);
    ulint count_sample(br_sample const& sample);
    ulint count_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples);
    std::vector<ulint> locate_sample(br_sample const& sample);
    std::vector<ulint> locate_samples(std::vector<br_sample> const& samples);
    std::vector<ulint> locate_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples);

    ulint count(std::string const& pattern); // count(P)
    std::vector<ulint> locate(std::string const& pattern); // locate(P)
    ulint maximal_exact_match(std::string const& pattern); // MEMs
    ulint full_task(ulint k, ulint t); // #substrings with length<=k & frequency>=t

    // suffix tree operations
    inline br_sample root() { return get_initial_sample(); }
    br_sample parent(br_sample const& sample);
    br_sample child(br_sample const& sample, uchar c);
    br_sample slink(br_sample const& sample);
    br_sample wlink(br_sample const& sample, uchar c);
    br_sample lca(br_sample const& left, br_sample const& right);
    ulint sdepth(br_sample const& sample);
    br_sample fchild(br_sample const& sample);
    br_sample nsibling(br_sample const& sample);
    std::vector<br_sample> children(br_sample const& sample);
    std::vector<uchar> child_chars(br_sample const& sample);
    uchar letter(br_sample const& sample, ulint i);
    std::string label(br_sample const& sample);
    bool ancestor(br_sample const& v, br_sample const& w) const { return v.contains(w); }
    inline bool is_leaf(br_sample const& sample) const { return sample.is_leaf(); }


    // I/O
    ulint serialize(std::ostream& out);
    void load(std::istream& in);
    void load(std::istream& in, ulint bl);
    void save_to_file(std::string const& path_prefix);
    void load_from_file(std::string const& path);
    ulint print_space(ulint fix);
    ulint print_space() { return print_space(length); }
    ulint get_space();


    // stats
    inline ulint text_size() { return bwt.size()-1; }
    inline ulint bwt_size() { return bwt.size(); }
    inline ulint border_length() { return length; }
    ulint number_of_runs(bool reversed=false) { return reversed ? bwtR.number_of_runs() : bwt.number_of_runs(); } // get the number of runs in the BWT
    uchar bwt_at(ulint i, bool reversed=false) { return reversed ? remap_inv[bwtR[i]] : remap_inv[bwt[i]]; }
    std::string get_bwt(bool reversed=false);
    


private:
    std::tuple<std::string, std::vector<range_t>, std::vector<range_t> > 
    sufsort(sdsl::int_vector<8>& text, sdsl::int_vector_buffer<>& sa);

    ulint _full_task_dfs(ulint k, ulint t, br_sample const& sample);

    static const uchar TERMINATOR = 1;
    
    // alphabet remapper
    std::vector<uchar> remap;
    std::vector<uchar> remap_inv;
    ulint sigma;

    // accumulated number of characters in lex order
    std::vector<ulint> F;
    
    // RLBWT
    rle_string_t bwt;
    ulint last_SA_val = 0;
    ulint r = 0;

    // RLBWT^R
    rle_string_t bwtR;
    ulint rR = 0;

    // left_extension
    sdsl::int_vector<> samples_first;
    sdsl::int_vector<> samples_last;
    
    // Phi (SA[i] -> SA[i-1])
    sparse_bitvector_t first;
    sdsl::int_vector<> first_to_run;
    
    // PhiI (SA[i] -> SA[i+1])
    sparse_bitvector_t last;
    sdsl::int_vector<> last_to_run;

    // right_extension
    sdsl::int_vector<> samples_firstR;
    sdsl::int_vector<> samples_lastR;

    // PhiR
    sparse_bitvector_t firstR;
    sdsl::int_vector<> first_to_runR;

    // PhiIR
    sparse_bitvector_t lastR;
    sdsl::int_vector<> last_to_runR;

    // determining the end of locate & left_contraction
    permuted_lcp<> plcp;

    // right_contraction
    permuted_lcp<> plcpR;

    // fixed length (parameter bl)
    ulint length;

    // left_contraction shorcut for |P|<=bl+1
    std::vector<sparse_bitvector_t> kmer;

    // right_contraction shortcut for |P|<=bl+1
    std::vector<sparse_bitvector_t> kmerR;

};

}; // namespace bri

#endif // INCLUDED_BR_INDEX_FULL_HPP
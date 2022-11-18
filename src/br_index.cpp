/*
 * fully functional bidirectional r-index
 * author: Yuma Arakawa
 */

#include "br_index.hpp"

namespace bri {

template<class T>
void print_vec(std::vector<T> const& vec) {
    for (auto x : vec) {
        std::cout << (ulint) x << " ";
    } std::cout << std::endl;
}

/*
* constructor. 
* \param input: string on which br-index is built
* \param length: parameter bl which effect the contraction algorithms and the index size.
* \param sais: flag determining if we use SAIS for suffix sort. 
*              otherwise we use divsufsort
* not suitable for texts of some gigabytes.
* use Prefix-Free Parsing builder instead. (pfpbuild.cpp)
*/
br_index::br_index(std::string const& input, ulint length, bool sais) {
    assert(length > 0);
    this->length = length;

    if (input.size() < 1)
    {

        std::cout << "Error: input string is empty" << std::endl;
        exit(1);

    }

    std::cout << "Text length = " << input.size() << std::endl;

    std::cout << "Parameter bl = " << length << std::endl << std::endl;

    std::cout << "(1/5) Remapping alphabet ... " << std::flush;
    
    // build RLBWT

    // configure & build indexes for sufsort & plcp
    sdsl::cache_config cc;

    // remap alphabet
    remap = std::vector<uchar>(256,0);
    remap_inv = std::vector<uchar>(256,0);
    {
        sigma = 1;
        std::vector<ulint> freqs(256,0);
        for (size_t i = 0; i < input.size(); ++i)
        {
            if (freqs[(uchar)input[i]]++ == 0) sigma++;
            if (sigma >= 255)
            {
                std::cout << "Error: alphabet cannot be remapped (overflow)" << std::endl;
                exit(1);
            }
        }
        uchar new_c = 2; // avoid reserved chars
        for (ulint c = 2; c < 256; ++c)
        {
            if (freqs[(uchar)c] != 0)
            {
                remap[(uchar)c] = new_c;
                remap_inv[new_c++] = (uchar)c;
            }
        }
    }

    std::cout << "done." << std::endl;
    std::cout << "(2/5) Building BWT, BWT^R, PLCP, PLCP^R and storing SA samples";
    if (sais) std::cout << " (SA-SAIS) ... " << std::flush;
    else std::cout << " (DIVSUFSORT) ... " << std::flush;

    // remap input text
    sdsl::int_vector<8> text(input.size());
    for (size_t i = 0; i < input.size(); ++i)
        text[i] = remap[(uchar)input[i]];

    sdsl::append_zero_symbol(text);

    // cache text
    sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, cc);
    sdsl::construct_config::byte_algo_sa = sais ? sdsl::SE_SAIS : sdsl::LIBDIVSUFSORT;
    
    // cache SA
    sdsl::construct_sa<8>(cc);
    // cache ISA 
    sdsl::construct_isa(cc);
    // cache LCP
    sdsl::construct_lcp_kasai<8>(cc);

    
    sdsl::int_vector_buffer<> sa(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
    last_SA_val = sa[sa.size()-1];
    auto bwt_and_samples = sufsort(text,sa);

    plcp = permuted_lcp<>(cc);

    // remove cache of text and SA
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, cc));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));




    // configure & build reversed indexes for sufsort
    sdsl::cache_config ccR;

    sdsl::int_vector<8> textR(input.size());
    for (ulint i = 0; i < input.size(); ++i)
        textR[i] = remap[(uchar)input[input.size()-1-i]];

    sdsl::append_zero_symbol(textR);

    // cache textR
    sdsl::store_to_cache(textR, sdsl::conf::KEY_TEXT, ccR);
    sdsl::construct_config::byte_algo_sa = sais ? sdsl::SE_SAIS : sdsl::LIBDIVSUFSORT;
    
    // cache SAR
    sdsl::construct_sa<8>(ccR);
    // cache ISAR
    sdsl::construct_isa(ccR);
    // cache LCPR
    sdsl::construct_lcp_kasai<8>(ccR);

    sdsl::int_vector_buffer<> saR(sdsl::cache_file_name(sdsl::conf::KEY_SA, ccR));
    auto bwt_and_samplesR = sufsort(textR,saR);

    // plcpR is needed for full index
    plcpR = permuted_lcp<>(ccR);

    // remove cache of textR and SAR
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, ccR));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, ccR));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, ccR));





    std::string& bwt_s = std::get<0>(bwt_and_samples);
    std::vector<range_t>& samples_first_vec = std::get<1>(bwt_and_samples);
    std::vector<range_t>& samples_last_vec = std::get<2>(bwt_and_samples);

    std::string& bwt_sR = std::get<0>(bwt_and_samplesR);
    std::vector<range_t>& samples_first_vecR = std::get<1>(bwt_and_samplesR);
    std::vector<range_t>& samples_last_vecR = std::get<2>(bwt_and_samplesR);

    std::cout << "done.\n(3/5) Run-length encoding BWT ... " << std::flush;


    // run length compression on BWT and BWTR
    bwt = rle_string_t(bwt_s);
    bwtR = rle_string_t(bwt_sR);

    // build F column (common between text and textR)
    F = std::vector<ulint>(256,0);

    for (uchar c : bwt_s) 
        F[c]++;

    for (ulint i = 255; i > 0; --i) 
        F[i] = F[i-1];

    F[0] = 0;

    for(ulint i = 1; i < 256; ++i) 
        F[i] += F[i-1];

    assert(input.size() + 1 == bwt.size());

    std::cout << "done." << std::endl << std::endl;


    r = bwt.number_of_runs();
    rR = bwtR.number_of_runs();

    assert(samples_first_vec.size() == r);
    assert(samples_last_vec.size() == r);

    assert(samples_first_vecR.size() == rR);
    assert(samples_last_vecR.size() == rR);

    int log_r = bitsize(r);
    int log_rR = bitsize(rR);
    int log_n = bitsize(bwt.size());

    std::cout << "Number of BWT equal-letter runs: r = " << r << std::endl;
    std::cout << "Rate n/r = " << double(bwt.size())/r << std::endl;
    std::cout << "log2(r) = " << std::log2(double(r)) << std::endl;
    std::cout << "log2(n/r) = " << std::log2(double(bwt.size())/r) << std::endl;

    std::cout << "Number of BWT^R equal-letter runs: rR = " << rR << std::endl << std::endl;

    // Phi, Phi inverse is needed only in forward case
    std::cout << "(4/5) Building predecessor/successor on R_c, R_c^R and Phi/Phi^{-1}/Phi_R/Phi^{-1}_R function ..." << std::flush;

    
    samples_last = sdsl::int_vector<>(r,0,log_n);
    samples_first = sdsl::int_vector<>(r,0,log_n);
    
    samples_firstR = sdsl::int_vector<>(rR,0,log_n);
    samples_lastR = sdsl::int_vector<>(rR,0,log_n);

    for (ulint i = 0; i < r; ++i)
    {
        samples_last[i] = samples_last_vec[i].first;
        samples_first[i] = samples_first_vec[i].first;
    }
    for (ulint i = 0; i < rR; ++i)
    {
        samples_lastR[i] = samples_last_vecR[i].first;
        samples_firstR[i] = samples_first_vecR[i].first;
    }

    // sort samples of first positions in runs according to text position
    std::sort(samples_first_vec.begin(), samples_first_vec.end());
    // sort samples of last positions in runs according to text position
    std::sort(samples_last_vec.begin(), samples_last_vec.end());

    // sort samples of first positions in runs according to text position
    std::sort(samples_first_vecR.begin(), samples_first_vecR.end());
    // sort samples of last positions in runs according to text position
    std::sort(samples_last_vecR.begin(), samples_last_vecR.end());

    // build Elias-Fano predecessor
    {
        std::vector<bool> first_bv(bwt_s.size(),false);
        for (auto p: samples_first_vec)
        {
            assert(p.first < first_bv.size());
            first_bv[p.first] = true;
        }
        first = sparse_bitvector_t(first_bv);
    }
    {
        std::vector<bool> last_bv(bwt_s.size(),false);
        for (auto p: samples_last_vec)
        {
            assert(p.first < last_bv.size());
            last_bv[p.first] = true;
        }
        last = sparse_bitvector_t(last_bv);
    }

    assert(first.rank(first.size()) == r);
    assert(last.rank(last.size()) == r);

    first_to_run = sdsl::int_vector<>(r,0,log_r);

    last_to_run = sdsl::int_vector<>(r,0,log_r);

    // construct first_to_run
    for (ulint i = 0; i < samples_first_vec.size(); ++i)
    {
        first_to_run[i] = samples_first_vec[i].second;
    }

    // construct last_to_run
    for (ulint i = 0; i < samples_last_vec.size(); ++i)
    {
        last_to_run[i] = samples_last_vec[i].second;
    }


    // build Elias-Fano predecessor for reversed direction
    {
        std::vector<bool> first_bvR(bwt_sR.size(),false);
        for (auto p: samples_first_vecR)
        {
            assert(p.first < first_bvR.size());
            first_bvR[p.first] = true;
        }
        firstR = sparse_bitvector_t(first_bvR);
    }
    {
        std::vector<bool> last_bvR(bwt_sR.size(),false);
        for (auto p: samples_last_vecR)
        {
            assert(p.first < last_bvR.size());
            last_bvR[p.first] = true;
        }
        lastR = sparse_bitvector_t(last_bvR);
    }

    assert(firstR.rank(firstR.size()) == rR);
    assert(lastR.rank(lastR.size()) == rR);

    first_to_runR = sdsl::int_vector<>(rR,0,log_rR);

    last_to_runR = sdsl::int_vector<>(rR,0,log_rR);

    // construct first_to_runR
    for (ulint i = 0; i < samples_first_vecR.size(); ++i)
    {
        first_to_runR[i] = samples_first_vecR[i].second;
    }

    // construct last_to_runR
    for (ulint i = 0; i < samples_last_vecR.size(); ++i)
    {
        last_to_runR[i] = samples_last_vecR[i].second;
    }

    std::cout << " done.\n(5/5) Building bitvectors for contraction shortcut S_l, E_l, S_l^R, E_l^R (l<=bl) ..." << std::flush;

    // construct kmer_start, kmer_end
    kmer = {};
    sdsl::int_vector<> lcp;
    sdsl::load_from_file(lcp, sdsl::cache_file_name(sdsl::conf::KEY_LCP, cc));
    for (ulint k = 0; k < length; ++k) {
        std::vector<bool> kmer_vector(text.size()+1,false);
        kmer_vector[0] = true;
        kmer_vector[text.size()] = true;

        bool yet;
        uchar tmp;

        yet = true;
        tmp = bwt_s[bwt_s.size()-1];
        for (ulint i = bwt_s.size()-1; i > 0; --i)
        {
            if (bwt_s[i-1] != tmp) {
                yet = true;
                tmp = bwt_s[i-1];
            }
            if (yet && lcp[i] <= k) {
                yet = false;
                kmer_vector[i] = true;
            }
        }
        yet = true;
        tmp = bwt_s[0];
        for (ulint i = 1; i < bwt_s.size(); ++i) {
            if (bwt_s[i] != tmp) {
                yet = true;
                tmp = bwt_s[i];
            }
            if (yet && lcp[i] <= k) {
                yet = false;
                kmer_vector[i] = true;
            }
        }

        kmer.push_back(sparse_bitvector_t(kmer_vector));
    }

    // construct kmerR
    kmerR = {};
    sdsl::int_vector<> lcpR;
    sdsl::load_from_file(lcpR, sdsl::cache_file_name(sdsl::conf::KEY_LCP, ccR));
    for (ulint k = 0; k < length; ++k) {
        std::vector<bool> kmer_vector(text.size()+1,false);
        kmer_vector[0] = true;
        kmer_vector[text.size()] = true;

        bool yet;
        uchar tmp;

        yet = true;
        tmp = bwt_sR[bwt_sR.size()-1];
        for (ulint i = bwt_sR.size()-1; i > 0; --i)
        {
            if (bwt_sR[i-1] != tmp) {
                yet = true;
                tmp = bwt_sR[i-1];
            }
            if (yet && lcpR[i] <= k) {
                yet = false;
                kmer_vector[i] = true;
            }
        }
        yet = true;
        tmp = bwt_sR[0];
        for (ulint i = 1; i < bwt_sR.size(); ++i) {
            if (bwt_sR[i] != tmp) {
                yet = true;
                tmp = bwt_sR[i];
            }
            if (yet && lcpR[i] <= k) {
                yet = false;
                kmer_vector[i] = true;
            }
        }

        kmerR.push_back(sparse_bitvector_t(kmer_vector));
    }

    // release LCP cache
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_LCP, cc));
    sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_LCP, ccR));

    std::cout << " done. " << std::endl << std::endl;
}

/*
 * get full BWT range
 */
range_t br_index::full_range()
{
    return {0,bwt_size()-1};
}

ulint br_index::LF(ulint i)
{
    auto c = bwt[i];
    return F[c] + bwt.rank(i,c);
}

/*
 * rn: BWT range of a string P
 * c:  remapped character
 * returns: BWT range of cP
 */
range_t br_index::LF(range_t rn, uchar c)
{

    if ((c == 255 && F[c] == bwt.size()) || F[c] >= F[c+1]) return {1,0};

    ulint c_before = bwt.rank(rn.first, c);

    ulint c_inside = bwt.rank(rn.second+1,c) - c_before;

    if (c_inside == 0) return {1,0};

    ulint lb = F[c] + c_before;

    return {lb, lb + c_inside - 1};

}

ulint br_index::LFR(ulint i)
{
    auto c = bwtR[i];
    return F[c] + bwtR.rank(i,c);
}

/*
 * rn: BWT^R range of a string P
 * c:  remapped character
 * returns: BWT^R range of cP
 */
range_t br_index::LFR(range_t rn, uchar c)
{

    if ((c == 255 && F[c] == bwt.size()) || F[c] >= F[c+1]) return {1,0};

    ulint c_before = bwtR.rank(rn.first, c);

    ulint c_inside = bwtR.rank(rn.second+1,c) - c_before;

    if (c_inside == 0) return {1,0};

    ulint lb = F[c] + c_before;

    return {lb, lb + c_inside - 1};

}

/*
 * inverse of LF (known as Psi)
 */
ulint br_index::FL(ulint i)
{

    // i-th character in first BWT column F
    auto c = F_at(i);

    // j: occurrences of c before i
    ulint j = i - F[c];

    return bwt.select(j,(uchar)c);

}

ulint br_index::FLR(ulint i)
{

    // i-th character in first BWT column F
    auto c = F_at(i);

    // j: occurrences of c before i
    ulint j = i - F[c];

    return bwtR.select(j,(uchar)c);
    
}

/*
 * Phi function
 * get SA[i] from SA[i+1]
 */
ulint br_index::Phi(ulint i)
{
    assert(i != bwt.size() - 1);

    ulint jr = first.predecessor_rank_circular(i);

    assert(jr <= r - 1);

    ulint k = first.select(jr);

    assert(jr < r - 1 || k == bwt.size() - 1);

    // distance from predecessor
    ulint delta = k < i ? i - k : i + 1;

    // check if Phi(SA[0]) is not called
    assert(first_to_run[jr] > 0);

    ulint prev_sample = samples_last[first_to_run[jr]-1];

    return (prev_sample + delta) % bwt.size();
}
/*
 * Phi inverse
 * get SA[i] from SA[i-1]
 */
ulint br_index::PhiI(ulint i)
{
    assert(i != last_SA_val);

    ulint jr = last.predecessor_rank_circular(i);

    assert(jr <= r - 1);

    ulint k = last.select(jr);

    assert(jr < r - 1 || k == bwt.size() - 1);

    // distance from predecessor
    ulint delta = k < i ? i - k : i + 1;

    // check if Phi(SA[0]) is not called
    assert(last_to_run[jr] < r-1);

    ulint prev_sample = samples_first[last_to_run[jr]+1];

    return (prev_sample + delta) % bwt.size();
}

/*
 * PhiR function
 * get SAR[i] from SAR[i+1]
 */
ulint br_index::PhiR(ulint i)
{
    assert(i != bwtR.size() - 1);

    ulint jr = firstR.predecessor_rank_circular(i);

    assert(jr <= rR - 1);

    ulint k = firstR.select(jr);

    assert(jr < rR - 1 || k == bwtR.size() - 1);

    // distance from predecessor
    ulint delta = k < i ? i - k : i + 1;

    // check if Phi(SA[0]) is not called
    assert(first_to_runR[jr] > 0);

    ulint prev_sample = samples_lastR[first_to_runR[jr]-1];

    return (prev_sample + delta) % bwtR.size();
}
/*
 * PhiR inverse
 * get SAR[i] from SAR[i-1]
 */
ulint br_index::PhiIR(ulint i)
{
    //assert(i != last_SA_val);

    ulint jr = lastR.predecessor_rank_circular(i);

    assert(jr <= rR - 1);

    ulint k = lastR.select(jr);

    assert(jr < rR - 1 || k == bwtR.size() - 1);

    // distance from predecessor
    ulint delta = k < i ? i - k : i + 1;

    // check if Phi(SA[0]) is not called
    assert(last_to_runR[jr] < rR-1);

    ulint prev_sample = samples_firstR[last_to_runR[jr]+1];

    return (prev_sample + delta) % bwtR.size();
}

/*
 * character of position i in column F
 */
uchar br_index::F_at(ulint i)
{

    ulint c = (std::upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
    assert(c < 256);
    assert(i >= F[c]);

    return (uchar)c;

}

/*
 * return BWT range of original char c (not remapped)
 */
range_t br_index::get_char_range(uchar c)
{
    // replace c with internal representation
    c = remap[c];

    if ((c == 255 && F[c] == bwt_size()) || F[c] >= F[c+1]) return {1,0};

    ulint lb = F[c];
    ulint rb = bwt_size() - 1;

    if (c < 255) rb = F[c+1] - 1;

    return {lb,rb};

}

/*
 * get a sample corresponding to an empty string
 */
br_sample br_index::get_initial_sample()
{
    return br_sample(full_range(), // entire SA range
                    full_range(), // entire SAR range
                    (samples_last[r-1]+1) % bwt.size(), // arbitrary sample
                    0,            // offset 0
                    0);           // null pattern

}

inline br_sample br_index::get_invalid_sample()
{
        return br_sample({1,0},   // empty SA range
                    {1,0},        // empty SAR range
                    0,            // arbitrary sample
                    0,            // offset 0
                    0);           // null pattern
}

/*
 * search the pattern cP (P:the current pattern)
 * returns SA&SAR range corresponding to cP
 * 
 * full version, update ranges for SA and SAR simultaneously
 * 
 * assumes c is original char (not remapped)
 */
br_sample br_index::left_extension(br_sample const& prev_sample, uchar c)
{
    // replace c with internal representation
    c = remap[c];

    br_sample sample(prev_sample);

    // get SA range of cP
    sample.range = LF(prev_sample.range,c);

    // pattern cP was not found
    if (sample.is_invalid()) return sample;

    // cP and aP occurs for some a s.t. a != c
    if (prev_sample.range.second - prev_sample.range.first != 
        sample.range.second      - sample.range.first)
    {
        
        // accumulated occ of aP (for any a s.t. a < c)
        ulint acc = 0;
        for (ulint a = 1; a < c; ++a)
        {
            range_t smaller_range = LF(prev_sample.range,(uchar)a);
            acc += (smaller_range.second+1) - smaller_range.first;
        }
        // get SAR range of (cP)^R
        sample.rangeR.second = sample.rangeR.first + acc + sample.range.second - sample.range.first;
        sample.rangeR.first = sample.rangeR.first + acc;


        // fint last c in range and get its sample
        // there must be at least one c due to the previous if clause
        ulint rnk = bwt.rank(prev_sample.range.second+1,c);
        assert(rnk > 0);

        // update p by corresponding BWT position
        ulint p = bwt.select(rnk-1,c);
        assert(p >= prev_sample.range.first && p <= prev_sample.range.second);

        // run number of position p
        ulint run_of_p = bwt.run_of_position(p);

        // update j by SA[p]
        if (bwt[prev_sample.range.second] == c)
            sample.j = samples_first[run_of_p];
        else
            sample.j = samples_last[run_of_p];

        // reset d
        sample.d = 0;

    }
    else // only c precedes P 
    {
        // increment offset by 1
        sample.d++;
    }
    sample.len++;
    return sample;
}

/*
 * search the pattern Pc (P:the current pattern)
 * return SAR&SA range corresponding to Pc
 * 
 * full version, updates ranges for SAR and SA simultaneously
 * 
 * assumes c is original char (not remapped)
 */
br_sample br_index::right_extension(br_sample const& prev_sample, uchar c)
{
    // replace c with internal representation
    c = remap[c];

    br_sample sample(prev_sample);

    // get SAR range of Pc
    sample.rangeR = LFR(prev_sample.rangeR,c);

    // pattern Pc was not found
    if (sample.is_invalid()) return sample;

    // Pc and Pa occurs for some a s.t. a != c
    if (prev_sample.rangeR.second - prev_sample.rangeR.first != 
        sample.rangeR.second      - sample.rangeR.first)
    {

        // accumulated occ of Pa (for any a s.t. a < c)
        ulint acc = 0;
        for (ulint a = 1; a < c; ++a)
        {
            range_t smaller_rangeR = LFR(prev_sample.rangeR,(uchar)a);
            acc += (smaller_rangeR.second+1) - smaller_rangeR.first;
        }
        // get SA range of Pc
        sample.range.second = sample.range.first + acc + sample.rangeR.second - sample.rangeR.first; 
        sample.range.first = sample.range.first + acc;


        // fint last c in range and get its sample
        // there must be at least one c due to the previous if clause
        ulint rnk = bwtR.rank(prev_sample.rangeR.second+1,c);
        assert(rnk > 0);

        // update pR by corresponding BWTR position
        ulint pR = bwtR.select(rnk-1,c);
        assert(pR >= prev_sample.rangeR.first && pR <= prev_sample.rangeR.second);

        // run number of position pR
        ulint run_of_pR = bwtR.run_of_position(pR);

        // j = SA[p]
        if (bwtR[prev_sample.rangeR.second] == c)
            sample.j = bwt.size()-2-samples_firstR[run_of_pR];
        else
            sample.j = bwt.size()-2-samples_lastR[run_of_pR];

        // reset d
        sample.d = sample.len;
    }
    /*else // only c follows P
    {
        // increment offset by 1, unncecessary for locate
        sample.dR++; 
    }*/
    sample.len++;
    return sample;
}

/*
 * match the pattern P from the current pattern cP
 * return SAR&SA range corresponding to P
 */
br_sample br_index::left_contraction(br_sample const& prev_sample)
{
    br_sample sample(prev_sample);
    assert(sample.len >= 1);

    if (sample.len == 1) return get_initial_sample();

    // get Psi(s), Psi(e)
    ulint psi_s = FL(sample.range.first);
    ulint psi_e = FL(sample.range.second);

    uchar c = bwt[psi_s];


    // updating range, rangeR
    ulint run_psi_s = bwt.run_of_position(psi_s);
    ulint run_start = bwt.run_start(run_psi_s);
    ulint run_psi_e = bwt.run_of_position(psi_e);
    ulint run_end = bwt.run_end(run_psi_e);
    assert(run_start <= psi_s);
    assert(psi_e <= run_end);

    if (sample.len <= length)
    {
        // predecessor in [0,psi_s]
        // if (run_start == psi_s) sample.range.first = kmer_start[sample.len-1].select(kmer_start[sample.len-1].rank(psi_s+1) - 1);
        if (run_start == psi_s) sample.range.first = kmer[sample.len-2].select(kmer[sample.len-2].rank(psi_s+1) - 1);
        else sample.range.first = psi_s;
        
        // successor in [psi_e,n-1]
        // if (run_end == psi_e) sample.range.second = kmer_end[sample.len-1].select(kmer_end[sample.len-1].rank(psi_e));
        if (run_end == psi_e) sample.range.second = kmer[sample.len-2].select(kmer[sample.len-2].rank(psi_e+1)) - 1;
        else sample.range.second = psi_e;
    }
    else 
    {
        sample.range.first = psi_s;
        if (run_start == psi_s) {
            ulint pos = (samples_first[run_psi_s]+1) % bwt.size();

            while (sample.range.first > 0 && plcp[pos] >= sample.len-1)
            {
                sample.range.first--;
                pos = Phi(pos);
            }
        }
        sample.range.second = psi_e;
        if (run_end == psi_e) {
            ulint pos = (samples_last[run_psi_e]+1) % bwt.size();

            while (sample.range.second < bwt.size()-1)
            {
                pos = PhiI(pos);
                if (plcp[pos] < sample.len-1) break;
                sample.range.second++;
            }
        }
    }
    
    // accumulated occ of aP (for any a s.t. a < c)
    ulint acc = 0;
    for (ulint a = 1; a < c; ++a)
    {
        range_t smaller_range = LF(sample.range,(uchar)a);
        acc += (smaller_range.second+1) - smaller_range.first;
    }
    // get range for SAR
    sample.rangeR.first -= acc;
    sample.rangeR.second = sample.rangeR.first + sample.range.second - sample.range.first;


    // updating j, d, len (very simple for contraction)
    if (sample.d == 0) sample.j++;
    else sample.d--;

    sample.len--;

    return sample;
}
/*
 * match the pattern P from the current pattern Pc
 * return SAR&SA range corresponding to P
 */
br_sample br_index::right_contraction(br_sample const& prev_sample)
{
    br_sample sample(prev_sample);
    assert(sample.len >= 1);

    if (sample.len == 1) return get_initial_sample();

    // get PsiR(sR), Psi(eR)
    ulint psiR_sR = FLR(sample.rangeR.first);
    ulint psiR_eR = FLR(sample.rangeR.second);

    uchar c = bwtR[psiR_sR];

    // updating range, rangeR
    ulint run_psiR_sR = bwtR.run_of_position(psiR_sR);
    ulint run_start = bwtR.run_start(run_psiR_sR);
    ulint run_psiR_eR = bwtR.run_of_position(psiR_eR);
    ulint run_end = bwtR.run_end(run_psiR_eR);
    assert(run_start <= psiR_sR);
    assert(psiR_eR <= run_end);

    if (sample.len <= length)
    {
        // predecessor in [0,psiR_sR]
        //if (run_start == psiR_sR) sample.rangeR.first = kmer_startR[sample.len-1].select(kmer_startR[sample.len-1].rank(psiR_sR+1) - 1);
        if (run_start == psiR_sR) sample.rangeR.first = kmerR[sample.len-2].select(kmerR[sample.len-2].rank(psiR_sR+1)-1);
        else sample.rangeR.first = psiR_sR;
        
        // successor in [psiR_eR,n-1]
        //if (run_end == psiR_eR) sample.rangeR.second = kmer_endR[sample.len-1].select(kmer_endR[sample.len-1].rank(psiR_eR));
        if (run_end == psiR_eR) sample.rangeR.second = kmerR[sample.len-2].select(kmerR[sample.len-2].rank(psiR_eR+1)) - 1;
        else sample.rangeR.second = psiR_eR;
    }
    else
    {
        sample.rangeR.first = psiR_sR;
        if (run_start == psiR_sR) {
            ulint pos = (samples_firstR[run_psiR_sR]+1) % bwtR.size();

            while (sample.rangeR.first > 0 && plcpR[pos] >= sample.len-1)
            {
                sample.rangeR.first--;
                pos = PhiR(pos);
            }
        }
        sample.rangeR.second = psiR_eR;
        if (run_end == psiR_eR) {
            ulint pos = (samples_lastR[run_psiR_eR]+1) % bwtR.size();

            while (sample.rangeR.second < bwtR.size()-1)
            {
                pos = PhiIR(pos);
                if (plcpR[pos] < sample.len-1) break;
                sample.rangeR.second++;
            }
        }
    }

    // accumulated occ of Pa (for any a s.t. a < c)
    ulint acc = 0;
    for (ulint a = 1; a < c; ++a)
    {
        range_t smaller_rangeR = LFR(sample.rangeR,(uchar)a);
        acc += (smaller_rangeR.second+1) - smaller_rangeR.first;
    }
    // get range for SA
    sample.range.first -= acc;
    sample.range.second = sample.range.first + sample.rangeR.second - sample.rangeR.first;

    // updating j, d, len (very simple for contraction)
    if (sample.d == sample.len - 1) 
    {
        sample.j--; sample.d--;
    }
    sample.len--;

    return sample;
}

// search the pattern
br_sample br_index::search(std::string const& pattern)
{
    br_sample sample(get_initial_sample());
    for (ulint i = 0; i < pattern.size(); ++i)
    {
        sample = right_extension(sample, pattern[i]);
        if (sample.is_invalid()) return sample;
    }
    return sample;
}

/*
 * count occurrences of current pattern P
 */
ulint br_index::count_sample(br_sample const& sample)
{
    return (sample.range.second + 1) - sample.range.first;
}

ulint br_index::count_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples)
{
    ulint res = 0;
    for (auto it = samples.begin(); it != samples.end(); ++it)
    {
        res += count_sample(it->second);
    }
    return res;
}

/*
 * locate occurrences of current pattern P
 * return them as std::vector
 * (space consuming if result is big)
 */
std::vector<ulint> br_index::locate_sample(br_sample const& sample)
{
    assert(sample.j >= sample.d);

    ulint sa = sample.j - sample.d;
    ulint pos = sa;

    std::vector<ulint> res;
    res.reserve(sample.range.second + 1 - sample.range.first);

    res.push_back(pos);

    while (plcp[pos] >= sample.len) 
    {
        pos = Phi(pos);
        res.push_back(pos);
    }
    pos = sa;
    while (true)
    {
        if (pos == last_SA_val) break;
        pos = PhiI(pos);
        if (plcp[pos] < sample.len) break;
        res.push_back(pos);
    }

    return res;
}

std::vector<ulint> br_index::locate_samples(std::vector<br_sample> const& samples)
{
    std::vector<ulint> res;
    for (auto s: samples)
    {
        auto tmp = locate_sample(s);
        res.insert(res.end(),tmp.begin(),tmp.end());
    }
    return res;
}

std::vector<ulint> br_index::locate_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples)
{
    std::vector<ulint> res;
    for (auto it = samples.begin(); it != samples.end(); ++it)
    {
        auto tmp = locate_sample(it->second);
        res.insert(res.end(),tmp.begin(),tmp.end());
    }
    return res;
}

/*
 * count the number of a given pattern
 */
ulint br_index::count(std::string const& pattern)
{
    br_sample sample(get_initial_sample());
    for (size_t i = 0; i < pattern.size(); ++i)
    {
        sample = right_extension(sample, pattern[i]);
        if (sample.is_invalid()) return 0;
    }
    return count_sample(sample);
}

/*
 * locate occurrences of a given pattern
 */
std::vector<ulint> br_index::locate(std::string const& pattern)
{
    br_sample sample(get_initial_sample());
    for (size_t i = 0; i < pattern.size(); ++i)
    {
        sample = right_extension(sample, pattern[i]);
        if (sample.is_invalid()) return {};
    }
    return locate_sample(sample);
}

// gets MEMs (returns max length of maximal substrings)
ulint br_index::maximal_exact_match(std::string const& pattern)
{
    ulint m = pattern.size();
    ulint j = 0, l = 0, max_l = 0;
    bool extended = false;

    br_sample init(get_initial_sample());
    br_sample sample(init);

    for (ulint i = 0; i < m; ++i)
    {
        while (j < m)
        {
            br_sample new_sample = right_extension(sample, (uchar)pattern[j]);
            if (new_sample.is_invalid()) break;
            extended = true;
            sample = new_sample;
            l++;
            j++;
        }
        if (extended)
        {
            if (l > max_l) max_l = l;
        }
        
        if (i == j)
        {
            sample = init;
            j++;
            l = 0;
        }
        else 
        {
            sample = left_contraction(sample);
            l--;
        }
        extended = false;
    }
    return max_l;
}

ulint br_index::_full_task_dfs(ulint k, ulint t, br_sample const& prev)
{
    if (prev.len >= k) return 0;
    
    ulint total = 0, acc = 0;

    for (ulint a = 2; a <= sigma; ++a) {
        br_sample sample(prev);

        sample.rangeR = LFR(prev.rangeR,(uchar)a);
        sample.len++;

        ulint cnt = sample.rangeR.second + 1 - sample.rangeR.first;
        if (cnt >= t) 
        {
            total += cnt;
            total += _full_task_dfs(k,t,sample);
        }
    }

    return total;
}

// compute total number of substrings with length<=k & frequency>=t
ulint br_index::full_task(ulint k, ulint t) {
    br_sample sample(get_initial_sample());
    return _full_task_dfs(k,t,sample);
}

// suffix tree op: parent
br_sample br_index::parent(br_sample const& sample)
{
    assert(sample.len > 0);

    br_sample par_sample(right_contraction(sample));

    if (sample.is_leaf() || par_sample.len == 0) return par_sample;

    br_sample tmp(sample);
    while (par_sample.size() == tmp.size())
    {
        tmp = par_sample;
        par_sample = right_contraction(tmp);
    }
    return par_sample;
}

// suffix tree op: child
br_sample br_index::child(br_sample const& sample, uchar c)
{
    if (sample.is_leaf()) return get_invalid_sample();

    br_sample new_sample(right_extension(sample,c));
    if (new_sample.is_invalid() || new_sample.is_leaf()) return new_sample;

    uchar a = bwtR[new_sample.rangeR.first];
    range_t rangeR = LFR(new_sample.rangeR,a);
    while (rangeR.second + 1 - rangeR.first == new_sample.rangeR.second + 1 - new_sample.rangeR.first)
    {
        new_sample.rangeR.first = rangeR.first;
        new_sample.rangeR.second = rangeR.second;
        new_sample.len++;
        rangeR = LFR(rangeR, bwtR[rangeR.first]);
    }
    return new_sample;
}

// suffix tree op: suffix-link
br_sample br_index::slink(br_sample const& sample)
{
    assert(sample.size() > 1);
    return left_contraction(sample);
}

// suffix tree op: weiner-link
br_sample br_index::wlink(br_sample const& sample, uchar c)
{
    return left_extension(sample,c);
}

// suffix tree op: lowest common ancestor
br_sample br_index::lca(br_sample const& left, br_sample const& right)
{
    if (left.len <= right.len) 
    {
        br_sample anc(left);
        while (!anc.contains(right)) anc = parent(anc);

        return anc;
    }
    else 
    {
        br_sample anc(right);
        while (!anc.contains(left)) anc = parent(anc);

        return anc;
    }
}

// suffix tree op: string-depth(v)
ulint br_index::sdepth(br_sample const& sample)
{ 
    if (sample.is_leaf()) return bwt.size() - (sample.j - sample.d);
    return sample.len; 
}
    
// suffix tree op: first-child(v)
br_sample br_index::fchild(br_sample const& sample)
{
    for (ulint a = 1; a <= sigma; ++a)
    {
        br_sample new_sample(child(sample,remap_inv[a]));
        if (!new_sample.is_invalid()) return new_sample;
    }
    return get_invalid_sample();
}
    
// suffix tree op: next-sibling(v)
br_sample br_index::nsibling(br_sample const& sample)
{

    br_sample par_sample(parent(sample));
    ulint p = sample.rangeR.first;
    for (ulint i = 0; i < sample.len-par_sample.len; ++i) p = FLR(p);

    uchar c = bwtR[p];
    br_sample new_sample;
    for (ulint a = c+1; a <= sigma; ++a)
    {
        new_sample = child(par_sample,remap_inv[a]);
        if (!new_sample.is_invalid()) return new_sample;
    }
    return get_invalid_sample();
}


// suffix tree op: children(v)
std::vector<br_sample> br_index::children(br_sample const& sample)
{
    std::vector<br_sample> res;
    for (ulint a = 1; a <= sigma; ++a)
    {
        br_sample tmp(child(sample,remap_inv[a]));
        if (!tmp.is_invalid()) res.push_back(tmp);
    }
    return res;
}

// child characters
std::vector<uchar> br_index::child_chars(br_sample const& sample)
{
    std::vector<uchar> res;
    for (ulint a = 1; a <= sigma; ++a)
    {
        range_t rangeR = LFR(sample.rangeR,a);
        if (rangeR.second+1-rangeR.first > 0) res.push_back(remap_inv[a]);
    }
    return res;
}

// suffix tree op: letter(v,i)
uchar br_index::letter(br_sample const& sample, ulint i)
{
    assert(i < sample.len || (sample.is_leaf() && i < bwt.size() - (sample.j-sample.d)));
    if ((i<<1) <= sample.len)
    {
        ulint p = sample.range.first;
        for (ulint j = 0; j <= i; ++j)
        {
            p = FL(p);
        }
        return remap_inv[bwt[p]];
    }
    else if (i < sample.len)
    {
        ulint p = sample.rangeR.first;
        for (ulint j = 0; j < sample.len-i; ++j)
        {
            p = FLR(p);
        }
        return remap_inv[bwtR[p]];
    }
    else // sample.is_leaf() && i < bwt.size() - (sample.j-sample.d)
    {
        ulint p = sample.rangeR.first;
        for (ulint k = 0; k < i - sample.len; ++k)
        {
            p = LFR(p);
        }
        return remap_inv[bwtR[p]];
    }
}

// suffix tree ops not supported:
// - tree-depth
// - level-ancestor-string
// - level-ancestor-tree
// they can be naively computed by iteratively applying right-contraction


// node label
std::string br_index::label(br_sample const& sample)
{
    std::string res;
    ulint p = sample.range.first;
    for (ulint i = 0; i < sample.len; ++i)
    {
        p = FL(p);
        res.push_back(remap_inv[(uchar)bwt[p]]);
    }
    return res;
}


ulint br_index::serialize(std::ostream& out)
{
    ulint w_bytes = 0;

    out.write((char*)&sigma,sizeof(sigma));
    out.write((char*)&length,sizeof(length));

    out.write((char*)remap.data(),256*sizeof(uchar));
    out.write((char*)remap_inv.data(),256*sizeof(uchar));

    out.write((char*)&last_SA_val,sizeof(last_SA_val));
    out.write((char*)F.data(),256*sizeof(ulint));

    w_bytes += sizeof(sigma)
                + sizeof(length)
                + 256*sizeof(uchar)
                + 256*sizeof(uchar)
                + sizeof(last_SA_val)
                + 256*sizeof(ulint);
    
    w_bytes += bwt.serialize(out);
    w_bytes += bwtR.serialize(out);

    w_bytes += samples_first.serialize(out);
    w_bytes += samples_last.serialize(out);

    w_bytes += first.serialize(out);
    w_bytes += first_to_run.serialize(out);

    w_bytes += last.serialize(out);
    w_bytes += last_to_run.serialize(out);

    w_bytes += samples_firstR.serialize(out);
    w_bytes += samples_lastR.serialize(out);

    w_bytes += firstR.serialize(out);
    w_bytes += first_to_runR.serialize(out);

    w_bytes += lastR.serialize(out);
    w_bytes += last_to_runR.serialize(out);

    w_bytes += plcp.serialize(out);
    w_bytes += plcpR.serialize(out);

    for (ulint k = 0; k < length; ++k)
    {
        w_bytes += kmer[k].serialize(out);

        w_bytes += kmerR[k].serialize(out);
    }


    return w_bytes;

}

void br_index::load(std::istream& in)
{

    in.read((char*)&sigma,sizeof(sigma));
    in.read((char*)&length,sizeof(length));

    remap = std::vector<uchar>(256);
    in.read((char*)remap.data(),256*sizeof(uchar));
    remap_inv = std::vector<uchar>(256);
    in.read((char*)remap_inv.data(),256*sizeof(uchar));
    
    in.read((char*)&last_SA_val,sizeof(last_SA_val));
    
    F = std::vector<ulint>(256);
    in.read((char*)F.data(),256*sizeof(ulint));

    bwt.load(in);
    bwtR.load(in);
    r = bwt.number_of_runs();
    rR = bwtR.number_of_runs();

    samples_first.load(in);
    samples_last.load(in);

    first.load(in);
    first_to_run.load(in);

    last.load(in);
    last_to_run.load(in);

    samples_firstR.load(in);
    samples_lastR.load(in);

    firstR.load(in);
    first_to_runR.load(in);

    lastR.load(in);
    last_to_runR.load(in);
    
    plcp.load(in);
    plcpR.load(in);

    kmer.resize(length);
    kmerR.resize(length);
    for (ulint k = 0; k < length; ++k)
    {
        kmer[k].load(in);

        kmerR[k].load(in);
    }

}
void br_index::load(std::istream& in, ulint bl)
{

    in.read((char*)&sigma,sizeof(sigma));
    in.read((char*)&length,sizeof(length));

    if (bl > length) {
        perror("br_index load error. Given bl is larger than idx file's parameter.");
        exit(1);
    }

    length = bl;

    remap = std::vector<uchar>(256);
    in.read((char*)remap.data(),256*sizeof(uchar));
    remap_inv = std::vector<uchar>(256);
    in.read((char*)remap_inv.data(),256*sizeof(uchar));
    
    in.read((char*)&last_SA_val,sizeof(last_SA_val));
    
    F = std::vector<ulint>(256);
    in.read((char*)F.data(),256*sizeof(ulint));

    bwt.load(in);
    bwtR.load(in);
    r = bwt.number_of_runs();
    rR = bwtR.number_of_runs();

    samples_first.load(in);
    samples_last.load(in);

    first.load(in);
    first_to_run.load(in);

    last.load(in);
    last_to_run.load(in);

    samples_firstR.load(in);
    samples_lastR.load(in);

    firstR.load(in);
    first_to_runR.load(in);

    lastR.load(in);
    last_to_runR.load(in);
    
    plcp.load(in);
    plcpR.load(in);

    kmer.resize(length);
    kmerR.resize(length);
    for (ulint k = 0; k < length; ++k)
    {
        kmer[k].load(in);
        kmerR[k].load(in);
    }
}

/*
 * save index to "{path_prefix}.bri" file
 */
void br_index::save_to_file(std::string const& path_prefix)
{

    std::string path = path_prefix + ".bri";
    
    std::ofstream out(path);
    serialize(out);
    out.close();

}

/*
 * load index file from path
 */
void br_index::load_from_file(std::string const& path)
{

    std::ifstream in(path);
    load(in);
    in.close();

}

/*
 * get stats
 */
ulint br_index::print_space(ulint fix)
{

    std::cout << "text length           : " << bwt.size() << std::endl;
    std::cout << "alphabet size         : " << sigma << std::endl;
    std::cout << "number of runs in bwt : " << bwt.number_of_runs() << std::endl;
    std::cout << "numbef of runs in bwtR: " << bwtR.number_of_runs() << std::endl << std::endl;
    
    ulint tot_bytes = sizeof(sigma)
                    + sizeof(length)
                    + 256*sizeof(uchar)
                    + 256*sizeof(uchar)
                    + sizeof(last_SA_val)
                    + 256*sizeof(ulint);
    
    std::cout << "(parameter bl=" << length << ")" << std::endl << std::endl;
    
    tot_bytes += bwt.print_space();
    tot_bytes += bwtR.print_space();
    std::cout << "total space for BWT: " << tot_bytes << " bytes" << std::endl << std::endl;

    tot_bytes += plcp.print_space();
    tot_bytes += plcpR.print_space();

    std::ofstream out("/dev/null");

    ulint bytes = 0;

    
    bytes =  samples_first.serialize(out);
    tot_bytes += bytes;
    std::cout << "samples_first:  " << bytes << " bytes" << std::endl;

    bytes =  samples_last.serialize(out);
    tot_bytes += bytes;
    std::cout << "samples_last:   " << bytes << " bytes" << std::endl;


    bytes =  first.serialize(out);
    tot_bytes += bytes;
    std::cout << "first:          " << bytes << " bytes" << std::endl;

    bytes =  first_to_run.serialize(out);
    tot_bytes += bytes;
    std::cout << "first_to_run:   " << bytes << " bytes" << std::endl;


    bytes =  last.serialize(out);
    tot_bytes += bytes;
    std::cout << "last:           " << bytes << " bytes" << std::endl;

    bytes =  last_to_run.serialize(out);
    tot_bytes += bytes;
    std::cout << "last_to_run:    " << bytes << " bytes" << std::endl;


    bytes =  samples_firstR.serialize(out);
    tot_bytes += bytes;
    std::cout << "samples_firstR: " << bytes << " bytes" << std::endl;

    bytes =  samples_lastR.serialize(out);
    tot_bytes += bytes;
    std::cout << "samples_lastR:  " << bytes << " bytes" << std::endl;


    bytes =  firstR.serialize(out);
    tot_bytes += bytes;
    std::cout << "firstR:         " << bytes << " bytes" << std::endl;

    bytes =  first_to_runR.serialize(out);
    tot_bytes += bytes;
    std::cout << "first_to_runR:  " << bytes << " bytes" << std::endl;


    bytes =  lastR.serialize(out);
    tot_bytes += bytes;
    std::cout << "lastR:          " << bytes << " bytes" << std::endl;

    bytes =  last_to_runR.serialize(out);
    tot_bytes += bytes;
    std::cout << "last_to_runR:   " << bytes << " bytes" << std::endl;

    std::cout << "kmer, kmerR: ";
    ulint kmer_bytes = 0;
    for (ulint k = 0; k < fix; ++k)
    {
        bytes =  kmer[k].serialize(out);
        kmer_bytes += bytes;
        tot_bytes += bytes;

        bytes =  kmerR[k].serialize(out);
        kmer_bytes += bytes;
        tot_bytes += bytes;
    }
    std::cout << kmer_bytes << " bytes" << std::endl;


    std::cout << "<total space of br-index>: " << tot_bytes << " bytes" << std::endl << std::endl;
    std::cout << "<bits/symbol>            : " << (double) tot_bytes * 8 / (double) bwt.size() << std::endl;

    return tot_bytes;

}

/*
 * get space complexity
 */
ulint br_index::get_space()
{

    ulint tot_bytes = sizeof(sigma)
                    + sizeof(length)
                    + 256*sizeof(uchar)
                    + 256*sizeof(uchar)
                    + sizeof(last_SA_val)
                    + 256*sizeof(ulint);

    tot_bytes += bwt.get_space();
    tot_bytes += bwtR.get_space();

    tot_bytes += plcp.get_space();
    tot_bytes += plcpR.get_space();

    std::ofstream out("/dev/null");

    tot_bytes += samples_first.serialize(out);
    tot_bytes += samples_last.serialize(out);

    tot_bytes += first.serialize(out);
    tot_bytes += first_to_run.serialize(out);

    tot_bytes += last.serialize(out);
    tot_bytes += last_to_run.serialize(out);

    tot_bytes += samples_firstR.serialize(out);
    tot_bytes += samples_lastR.serialize(out);

    tot_bytes += firstR.serialize(out);
    tot_bytes += first_to_runR.serialize(out);

    tot_bytes += lastR.serialize(out);
    tot_bytes += last_to_runR.serialize(out);


    for (ulint k = 0; k < length; ++k)
    {
        tot_bytes += kmer[k].serialize(out);

        tot_bytes += kmerR[k].serialize(out);
    }

    return tot_bytes;

}

/*
 * get string representation of BWT
 */
std::string br_index::get_bwt(bool reversed)
{
    if (!reversed)
    {
        std::string res(bwt.to_string());
        for (size_t i = 0; i < res.size(); ++i)
            res[i] = remap_inv[(uchar)res[i]];
        return res;
    } else {
        std::string res(bwtR.to_string());
        for (size_t i = 0; i < res.size(); ++i)
            res[i] = remap_inv[(uchar)res[i]];
        return res;
    }
}

/*
 * builds BWT from SA
 */
std::tuple<std::string, std::vector<range_t>, std::vector<range_t> > 
br_index::sufsort(sdsl::int_vector<8>& text, sdsl::int_vector_buffer<>& sa)
{
    std::string bwt_s;
    std::vector<range_t> samples_first;
    std::vector<range_t> samples_last;

    {
        for (ulint i = 0; i < sa.size(); ++i)
        {
            auto x = sa[i];

            assert(x <= text.size());

            if (x > 0) 
                bwt_s.push_back((uchar)text[x-1]);
            else 
                bwt_s.push_back(TERMINATOR);
            
            // insert samples at beginnings of runs
            if (i > 0)
            {
                if (i==1 || (i>1 && bwt_s[i-1] != bwt_s[i-2]))
                {
                    samples_first.push_back({
                        sa[i-1] > 0
                        ? sa[i-1] - 1
                        : sa.size() - 1,
                        samples_first.size()
                    });
                }
                if (i==sa.size()-1 && bwt_s[i] != bwt_s[i-1])
                {
                    samples_first.push_back({
                        sa[i] > 0
                        ? sa[i] - 1
                        : sa.size() - 1,
                        samples_first.size()
                    });
                }
            }

            // insert samples at ends of runs
            if (i > 0)
            {
                if (bwt_s[i-1] != bwt_s[i])
                {
                    samples_last.push_back({
                        sa[i-1] > 0
                        ? sa[i-1] - 1
                        : sa.size() - 1,
                        samples_last.size()
                    });
                }
                if (i == sa.size()-1)
                {
                    samples_last.push_back({
                        sa[i] > 0
                        ? sa[i] - 1
                        : sa.size() - 1,
                        samples_last.size()
                    });
                }
            }
        }
    }

    return std::tuple<std::string, std::vector<range_t>, std::vector<range_t> >
        (bwt_s, samples_first, samples_last);
}

};
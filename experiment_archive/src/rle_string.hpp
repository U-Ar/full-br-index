/*
 * rle_string.hpp
 *
 *  Original author: nicola
 *
 *  A run-length encoded string with rank/access functionalities.
 *
 *
 *  space of the structure: R * (H0 + log(n/R) + log(n/R)/B ) (1+o(1)) bits, n being text length,
 *  R number of runs, B block length, and H0 zero-order entropy of the run heads.
 *
 *  Time for all operations: O( B*(log(n/R)+H0) )
 *
 *  From the paper
 *
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza and Mathieu Raffinot.
 *  Flexible Indexing of Repetitive Collections. Computability in Europe (CiE) 2017)
 *
 */

#ifndef INCLUDED_RLE_STRING_HPP
#define INCLUDED_RLE_STRING_HPP

#include "definitions.hpp"
#include "huffman_string.hpp"
#include "sparse_sd_vector.hpp"

namespace bri {

template<
    class sparse_bitvector_t = sparse_sd_vector,
    class string_t = huffman_string
>
class rle_string {

public:
    rle_string() {}

    /*
     * constructor
     * \param input the input string
     * \param B block size
     */
    rle_string(std::string& input, ulint B = 2)
    {
        assert(!contains0(input));

        this->B = B;
        n = input.size();
        r = 0;

        auto runs_per_letter_bv = std::vector<std::vector<bool> >(256);

        std::vector<bool> runs_bv;
        std::string run_heads_s;

        uchar last_c = input[0];

        for (ulint i = 1; i < input.size(); ++i)
        {
            if (uchar(input[i]) != last_c)
            {
                run_heads_s.push_back(last_c);
                runs_per_letter_bv[last_c].push_back(true);

                last_c = input[i];

                // push back a bit set only at the end of a block
                runs_bv.push_back(r%B == B-1);

                r++;

            } else {

                runs_bv.push_back(false);
                runs_per_letter_bv[last_c].push_back(false);

            }
        }

        run_heads_s.push_back(last_c);
        runs_per_letter_bv[last_c].push_back(true);
        runs_bv.push_back(false);
        r++;

        assert(run_heads_s.size()==r);
        assert(r==count_runs(input));
        assert(runs_bv.size()==input.size());

        ulint t = 0;
        for (ulint i = 0; i < 256; ++i) t += runs_per_letter_bv[i].size();
        assert(t == input.size());

        runs = sparse_bitvector_t(runs_bv);

        runs_per_letter = std::vector<sparse_bitvector_t>(256);
        for (ulint i = 0; i < 256; ++i)
            runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);

        run_heads = string_t(run_heads_s);

        assert(run_heads.size() == r);

    }

    uchar operator[](size_t i)
    {
        assert(i < n);
        return run_heads[run_of(i).first];
    }

    /*
     * index of i-th character c 
     */
    size_t select(ulint i, uchar c)
    {
        assert(i<runs_per_letter[c].size());

        // i-th c is inside j-th c-run
        ulint j = runs_per_letter[c].rank(i);

        assert(j==0 || i >= runs_per_letter[c].select(j-1)+1);

        // starting position of i-th c inside its run
        ulint before = (j==0 ? i : i-(runs_per_letter[c].select(j-1)+1));

        // position in run_heads
        ulint p = run_heads.select(j,c);

        // k = number of bits before position of interest in the main string
        // k is initialized looking at the sampled runs
        assert(p/B==0 || p/B - 1 < runs.number_of_1());
        ulint k = (p/B == 0 ? 0 : runs.select(p/B-1)+1);

        // add remaining run length to k
        for (ulint t = (p/B)*B; t < p; ++t) k += run_at(t);

        return k + before;
    }

    /*
     * number of c before position i
     */
    ulint rank(size_t i, uchar c)
    {

        assert(i <= n);

        // c does not exist
        if (runs_per_letter[c].size() == 0) return 0;

        // number of all c
        if (i == n) return runs_per_letter[c].size();

        ulint last_block = runs.rank(i);
        ulint current_run = last_block * B;

        // current position in the string
        ulint pos = 0;
        if (last_block > 0) pos = runs.select(last_block-1) + 1;

        assert(pos <= i);

        ulint dist = i - pos;

        // scan at most B runs
        while (pos < i) 
        {
            pos += run_at(current_run);
            current_run++;

            if (pos <= i) dist = i - pos;
        }

        if (pos > i) current_run--;

        assert(current_run < r);

        // number of c runs before the current run
        ulint rk = run_heads.rank(current_run, c);

        // number of c before i in the current run
        ulint tail = (run_heads[current_run] == c) * dist;

        if (rk == 0) return tail;

        return runs_per_letter[c].select(rk-1)+1+tail;

    }

    /*
     * run number of text position i
     */
    ulint run_of_position(size_t i)
    {

        assert(i < n);

        ulint last_block = runs.rank(i);
        ulint current_run = last_block * B;

        //current position in the string: the first of a block
        ulint pos = 0;
        if(last_block > 0)
            pos = runs.select(last_block-1)+1;

        assert(pos <= i);

        ulint dist = i - pos;

        //scan at most B runs
        while(pos < i){

            pos += run_at(current_run);
            current_run++;

            if(pos <= i) dist = i - pos;

        }

        if(pos>i) current_run--;

        assert(current_run < r);

        return current_run;

    }

    ulint size() { return n; }

    /*
     * inclusive range of j-th run
     */
    range_t run_range(ulint j)
    {

        assert(j < run_heads.size());

        ulint this_block = j/B;
        ulint current_run = this_block * B;
        ulint pos = (this_block == 0 ? 0 : runs.select(this_block-1)+1);

        while (current_run < j)
        {
            pos += run_at(current_run);
            current_run++;
        }

        assert(current_run == j);

        return {pos, pos + run_at(j) - 1};

    }

    ulint run_start(ulint j)
    {
        assert(j < run_heads.size());

        ulint this_block = j/B;
        ulint current_run = this_block * B;
        ulint pos = (this_block == 0 ? 0 : runs.select(this_block-1)+1);

        while (current_run < j)
        {
            pos += run_at(current_run);
            current_run++;
        }

        assert(current_run == j);

        return pos;
    }

    ulint run_end(ulint j)
    {
        assert(j < run_heads.size());

        ulint this_block = j/B;
        ulint current_run = this_block * B;
        ulint pos = (this_block == 0 ? 0 : runs.select(this_block-1)+1);

        while (current_run < j)
        {
            pos += run_at(current_run);
            current_run++;
        }

        assert(current_run == j);

        return pos + run_at(j) - 1;
    }


    /*
     * length of i-th run
     */
    ulint run_at(ulint i)
    {

        assert(i < r);
        uchar c = run_heads[i];

        return runs_per_letter[c].gap_at(run_heads.rank(i,c));

    }

    ulint number_of_runs() { return r; }

    ulint serialize(std::ostream& out)
    {

        ulint w_bytes = 0;

        out.write((char*)&n,sizeof(n));
        out.write((char*)&r,sizeof(r));
        out.write((char*)&B,sizeof(B));

        w_bytes += sizeof(n) + sizeof(r) + sizeof(B);

        if (n == 0) return w_bytes;

        w_bytes += runs.serialize(out);

        for (ulint i = 0; i < 256; ++i)
            w_bytes += runs_per_letter[i].serialize(out);

        w_bytes += run_heads.serialize(out);

        return w_bytes;

    }

    void load(std::istream& in)
    {

        in.read((char*)&n,sizeof(n));
        in.read((char*)&r,sizeof(r));
        in.read((char*)&B,sizeof(B));

        if (n == 0) return;

        runs.load(in);

        runs_per_letter = std::vector<sparse_bitvector_t>(256);

        for (ulint i = 0; i < 256; ++i)
            runs_per_letter[i].load(in);
        
        run_heads.load(in);

    }

    std::string to_string()
    {
        
        std::string s;

        for (size_t i = 0; i < size(); ++i)
            s.push_back(operator[](i));
        
        return s;

    }

    ulint print_space()
    {

        ulint tot_bytes = 0;

        std::cout << "space for run-length encoded string:" << std::endl;

        std::ofstream out("/dev/null");
        auto bytesize = runs.serialize(out);
        tot_bytes += bytesize;
        std::cout << "- bitvector runs: " << bytesize << " bytes" << std::endl;

        bytesize = 0;
        for (auto r: runs_per_letter) bytesize += r.serialize(out);
        tot_bytes += bytesize;
        std::cout << "- bitvectors runs-per-letter: " << bytesize << " bytes" << std::endl;

        bytesize = run_heads.serialize(out);
        tot_bytes += bytesize;
        std::cout << "- run heads: " << bytesize << " bytes" << std::endl;

        return tot_bytes;

    }

    ulint get_space()
    {
        ulint tot_bytes = 0;

        std::ofstream out("/dev/null");
        auto bytesize = runs.serialize(out);
        tot_bytes += bytesize;

        bytesize = 0;
        for (auto r: runs_per_letter) bytesize += r.serialize(out);
        tot_bytes += bytesize;

        bytesize = run_heads.serialize(out);
        tot_bytes += bytesize;

        return tot_bytes;

    }

    // <j(run number of position i), p(last position of j-th run)>
    std::pair<ulint,ulint> run_of(ulint i)
    {
        ulint last_block = runs.rank(i);
        ulint current_run = last_block * B;

        // current position in the string
        ulint pos = 0;
        if (last_block > 0) pos = runs.select(last_block-1) + 1;

        assert(pos <= i);

        while (pos < i)
        {
            pos += run_at(current_run);
            current_run++;
        }

        assert(pos >= i);

        if (pos > i)
        {
            current_run--;
        } else {
            pos += run_at(current_run);
        }

        assert(pos > 0);
        assert(current_run < r);

        return {current_run, pos-1};
    }

private:

    // static member func to count the number of runs in s
    static ulint count_runs(std::string& s)
    {
        ulint runs = 1;

        for (ulint i = 1; i < s.size(); ++i)
        {
            if (s[i] != s[i-1]) runs++;
        }

        return runs;
    }

    bool contains0(std::string& s)
    {
        for (auto c: s) if (c==0) return true;
        return false;
    }

    ulint B = 0;

    sparse_bitvector_t runs;

    // runs for each letter
    std::vector<sparse_bitvector_t> runs_per_letter;

    // run heads with rank/select support
    string_t run_heads;

    // text length
    ulint n = 0;

    // number of runs
    ulint r = 0;

};

typedef rle_string<sparse_sd_vector> rle_string_sd;

};


#endif
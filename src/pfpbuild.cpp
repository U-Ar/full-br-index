#include <unistd.h>

#include "definitions.hpp"
#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include "permuted_lcp.hpp"
#include "br_index.hpp"
extern "C" {
    #include "pfpbuilder/utils.h"
}

using namespace bri;

template<class T>
void print_vec(std::vector<T> const& vec) {
    for (auto x : vec) {
        std::cout << (ulint) x << " ";
    } std::cout << std::endl;
}

// index file extension
const std::string EXTIDX = "bri";

// struct containing command line parameters and other globals
struct Args {
    std::string input_file = "";
    std::string output_base = "";
    int bl = 8;
    int verbose=0;
    bool inmemory = false;
    FILE *tmp_parse_file, *last_file, *sa_file; 
};

namespace bri {

class br_index_builder {
    br_index idx;
public:
    // input: file basename of Prefix-Free Parsing (.bwt, .ssa, .esa, .rev.bwt, .rev.ssa, .rev.esa are necessary)
    // bl: parameter
    br_index_builder& build_from_pfp(std::string const& input, int bl) {

        std::cout << "Building br-index on " << input << std::endl;
        std::cout << "  using PFP parsing files." << std::endl;

        idx = br_index();
        idx.length = (ulint)bl;

        std::cout << "bl: " << bl << std::endl;

        // read .bwt
        std::ifstream fbwt(input + ".bwt", std::ifstream::ate);
        if (!fbwt.is_open()) {
            perror(__func__);
            throw new std::runtime_error("Cannot open file " + input + ".bwt");
        }
        long size = fbwt.tellg();
        fbwt.seekg(0);

        std::cout << "Remapping alphabet ... " << std::flush;

        // construct alphabet remapper (null character \0 is mapped to \1)
        idx.remap = std::vector<uchar>(256,0); idx.remap[0] = 1;
        idx.remap_inv = std::vector<uchar>(256,0); idx.remap_inv[1] = 0;
        idx.sigma = 0;
        idx.F = std::vector<ulint>(256,0);
        {
            std::vector<ulint> freqs(256,0);
            for (size_t i = 0; i < size; ++i) {
                if (freqs[(uchar)fbwt.get()]++ == 0) idx.sigma++;
                if (idx.sigma >= 255) {
                    die("Alphabet overflow");
                }
            }
            assert(freqs[0]==1);
            uchar new_c = 2;
            for (ulint c = 1; c < 256; ++c) {
                if (freqs[(uchar)c] != 0) {
                    idx.remap[(uchar)c] = new_c;
                    idx.remap_inv[new_c++] = (uchar)c;
                }
            }
            // build F column
            for (ulint i = 0; i < 256; ++i) idx.F[i] = freqs[idx.remap_inv[i]];
            for (ulint i = 255; i > 0; --i) idx.F[i] = idx.F[i-1];
            idx.F[0] = 0;
            for (ulint i = 1; i < 256; ++i) idx.F[i] += idx.F[i-1];
        }
        assert(idx.remap[0] == 1);
        assert(idx.remap_inv[1] == 0);

        std::cout << "done.\nBuilding RLBWT ... " << std::flush;

        // build RLBWT with remapper
        fbwt.seekg(0);
        idx.bwt.load_from_plain(fbwt,size,idx.remap);
        fbwt.close();
        assert(size==idx.bwt.size());

        idx.r = idx.bwt.number_of_runs();
        ulint r = idx.r;
        int log_r = bitsize(r);
        int log_n = bitsize(size);

        std::cout << "done.\nNumber of runs in BWT  r: " << r << std::endl;

        std::cout << "Reading SA samples at run boundaries ... " << std::flush;

        // read .ssa
        FILE* file_ssa = open_aux_file(input.c_str(),EXTSSA,"rb");

        idx.samples_first = sdsl::int_vector<>(r,0,log_n);
        idx.first_to_run = sdsl::int_vector<>(r,0,log_r);
        {
            auto pos_run_pairs = std::vector<std::pair<ulint,ulint>>(r);
            /*ulint* sa_val;
            for (ulint i = 0; i < r; ++i) {
                fread(sa_val,SABYTES,1,file_ssa);
                size_t s = fread(sa_val,SABYTES,1,file_ssa);
                if (s!=1) die(".ssa read failed");
                idx.samples_first[i] = *sa_val > 0 ? *sa_val-1 : size-1;
                pos_run_pairs[i] = {*sa_val > 0 ? *sa_val-1 : size-1,i};
            }*/
            ulint* buf = new ulint[2];
            for (ulint i = 0; i < r; ++i) {
                size_t s = fread(buf, SABYTES, 2, file_ssa);
                if (s!=2) die(".ssa read failed.");
                ulint sa_val = buf[1] > 0 ? buf[1] - 1 : size-1;
                idx.samples_first[i] = sa_val;
                pos_run_pairs[i] = {sa_val,i};
            }
            delete[] buf;

            std::sort(pos_run_pairs.begin(),pos_run_pairs.end());

            std::vector<ulint> first_pos;
            for (ulint i = 0; i < r; ++i) {
                first_pos.push_back(pos_run_pairs[i].first);
                idx.first_to_run[i] = pos_run_pairs[i].second;
            }
            idx.first = br_index::sparse_bitvector_t(first_pos.cbegin(),first_pos.cend());
        }
        fclose(file_ssa);

        // read.esa
        FILE* file_esa = open_aux_file(input.c_str(),EXTESA,"rb");

        idx.samples_last = sdsl::int_vector<>(r,0,log_n);
        idx.last_to_run = sdsl::int_vector<>(r,0,log_r);
        {
            auto pos_run_pairs = std::vector<std::pair<ulint,ulint>>(r);
            /* ulint* sa_val;
            for (ulint i = 0; i < r; ++i) {
                fread(sa_val,SABYTES,1,file_esa);
                size_t s = fread(sa_val,SABYTES,1,file_esa);
                if (s!=1) die(".esa read failed");
                idx.samples_last[i] = *sa_val > 0 ? *sa_val-1 : size-1;
                pos_run_pairs[i] = {*sa_val > 0 ? *sa_val-1 : size-1,i};
            }
            idx.last_SA_val = *sa_val;*/
            ulint* buf = new ulint[2];
            for (ulint i = 0; i < r; ++i) {
                size_t s = fread(buf, SABYTES, 2, file_esa);
                if (s!=2) die(".esa read failed.");
                ulint sa_val = buf[1] > 0 ? buf[1] - 1 : size-1;
                idx.samples_last[i] = sa_val;
                pos_run_pairs[i] = {sa_val,i};
            }
            idx.last_SA_val = buf[1];
            delete[] buf;


            std::sort(pos_run_pairs.begin(),pos_run_pairs.end());

            std::vector<ulint> last_pos;
            for (ulint i = 0; i < r; ++i) {
                last_pos.push_back(pos_run_pairs[i].first);
                idx.last_to_run[i] = pos_run_pairs[i].second;
            }
            idx.last = br_index::sparse_bitvector_t(last_pos.cbegin(),last_pos.cend());
        }
        fclose(file_esa);

        std::cout << "done.\nBuilding run-length compressed PLCP ... " << std::flush;

        // construct run-length encoded PLCP
        {
            std::vector<ulint> ones, zeros;

            ulint i = idx.FL(0); // ISA[0]
            assert(i>0);
            uchar c = idx.F_at(i); // T[0]
            ulint i0 = i-1;
            uchar c0 = idx.F_at(i-1); // T[\phi(0)]
            ulint l = 0, prev_l = 0;
            ulint acc0 = 0, acc1 = 0;
            ones.push_back(0);
            for (ulint j = 0; j < size-1; ++j) {
                while (c == c0) {
                    // compute FL(i) & FL(i0), compare BWT chars
                    l++;
                    ulint p = i - idx.F[c];
                    i = idx.bwt.select(p,c);
                    c = idx.F_at(i);
                    ulint p0 = i0 - idx.F[c0];
                    i0 = idx.bwt.select(p0,c0);
                    c0 = idx.F_at(i0);
                }
                if (l==0) {
                    ulint p = i - idx.F[c];
                    i = idx.bwt.select(p,c); // i = FL(i)
                    c = idx.F_at(i);
                    assert(i>0 || j == size-2);
                    i0 = i - 1;
                    c0 = idx.F_at(i0);
                }

                if (j == 0) {
                    if (l) {
                        acc0 += l - 1;
                        zeros.push_back(acc0);
                        acc1++;
                    } else {
                        acc1++;
                    }
                }
                else if (l + 1 - prev_l) {
                    ones.push_back(acc1);
                    acc0 += l + 1 - prev_l;
                    zeros.push_back(acc0);
                    acc1++;
                } else {
                    acc1++;
                }

                prev_l = l;
                if (l) l--;      
            }
            ones.push_back(acc1);

            idx.plcp = permuted_lcp<>(size,ones,zeros);
        }

        std::cout << "done.\nBuilding kmer[0,bl) ... " << std::flush;

        // 各runの境界につき、その前後を探索する形？
        // run_range(run_num)でrunの範囲が分かる。そんでそこからSAの値もわかるので、
        // そこ起点で広げたらいけるのでは？という気がしている
        
        // consruct kmer[0,bl) (contraction shortcut bitvector)
        {
            auto kmer_pos = std::vector<std::vector<ulint>>(bl);
            // first run head
            for (ulint i = 0; i < bl; ++i) kmer_pos[i].push_back(0);

            ulint prev_run_start = 0;
            ulint run_start, run_end;
            ulint cur_sa, min_lcp_l, min_lcp_r;

            for (ulint i = 1; i < r; ++i) {
                auto run_range = idx.bwt.run_range(i);
                run_start = run_range.first;
                run_end = run_range.second;

                ulint start_sa = idx.samples_first[i] < size-1 ? idx.samples_first[i]+1 : 0;
                min_lcp_l = min_lcp_r = idx.plcp[start_sa];
                if (min_lcp_l > bl) min_lcp_l = min_lcp_r = bl;
                
                for (ulint l = min_lcp_l; l < bl; ++l) kmer_pos[l].push_back(run_start);

                // backward scan from run head
                cur_sa = start_sa;
                for (ulint p = run_start-1; p > prev_run_start; --p) {
                    if (min_lcp_l==0) break;

                    cur_sa = idx.Phi(cur_sa);
                    ulint lcp = idx.plcp[cur_sa];
                    if (lcp < min_lcp_l) {
                        for (ulint l = lcp; l < min_lcp_l; ++l) {
                            if (kmer_pos[l].back() != p) kmer_pos[l].push_back(p);
                        } 
                        min_lcp_l = lcp;
                    }
                }
                // forward scan from run head
                cur_sa = start_sa;
                for (ulint p = run_start+1; p <= run_end; ++p) {
                    if (min_lcp_r==0) break;

                    cur_sa = idx.PhiI(cur_sa);
                    ulint lcp = idx.plcp[cur_sa];
                    if (lcp < min_lcp_r) {
                        for (ulint l = lcp; l < min_lcp_r; ++l) kmer_pos[l].push_back(p);
                        min_lcp_r = lcp;
                    }
                }

                prev_run_start = run_start;
            }
            // last run tail
            for (ulint i = 0; i < bl; ++i) kmer_pos[i].push_back(size);

            idx.kmer = std::vector<br_index::sparse_bitvector_t>(bl);
            for (ulint i = 0; i < bl; ++i) {
                idx.kmer[i] = br_index::sparse_bitvector_t(kmer_pos[i].cbegin(),kmer_pos[i].cend());
            }
        }
        
        std::cout << "done." << std::endl;
        std::cout << "Start building components for the reversed direction." << std::endl;

        // load from reversed files
        std::string input_rev = input + ".rev";

        // read .rev.bwt
        std::ifstream fbwt_rev(input_rev + ".bwt");
        if (!fbwt_rev.is_open()) {
            perror(__func__);
            throw new std::runtime_error("Cannot open file " + input_rev + ".bwt");
        }

        std::cout << "Building RLBWT^R ... " << std::flush;

        // build RLBWT^R with remapper
        fbwt_rev.seekg(0);
        idx.bwtR.load_from_plain(fbwt_rev,size,idx.remap);
        fbwt_rev.close();
        assert(size==idx.bwtR.size());

        idx.rR = idx.bwtR.number_of_runs();
        ulint rR = idx.rR;
        int log_rR = bitsize(rR);

        std::cout << "done.\nNumber of runs in BWT^R  rR: " << rR << std::endl;

        std::cout << "Reading SA^R samples at run boundaries ... " << std::flush;

        // read .rev.ssa
        FILE* file_ssa_rev = open_aux_file(input_rev.c_str(),EXTSSA,"rb");

        idx.samples_firstR = sdsl::int_vector<>(rR,0,log_n);
        idx.first_to_runR = sdsl::int_vector<>(rR,0,log_rR);
        {
            auto pos_run_pairs = std::vector<std::pair<ulint,ulint>>(rR);
            /*ulint* sa_val;
            for (ulint i = 0; i < rR; ++i) {
                fread(sa_val,SABYTES,1,file_ssa_rev);
                size_t s = fread(sa_val,SABYTES,1,file_ssa_rev);
                if (s!=1) die(".rev.ssa read failed");
                idx.samples_firstR[i] = *sa_val > 0 ? *sa_val-1 : size-1;
                pos_run_pairs[i] = {*sa_val > 0 ? *sa_val-1 : size-1,i};
            }*/
            ulint* buf = new ulint[2];
            for (ulint i = 0; i < rR; ++i) {
                size_t s = fread(buf, SABYTES, 2, file_ssa_rev);
                if (s!=2) die(".rev.ssa read failed.");
                ulint sa_val = buf[1] > 0 ? buf[1] - 1 : size-1;
                idx.samples_firstR[i] = sa_val;
                pos_run_pairs[i] = {sa_val,i};
            }
            delete[] buf;

            std::sort(pos_run_pairs.begin(),pos_run_pairs.end());

            std::vector<ulint> first_pos;
            for (ulint i = 0; i < rR; ++i) {
                first_pos.push_back(pos_run_pairs[i].first);
                idx.first_to_runR[i] = pos_run_pairs[i].second;
            }
            idx.firstR = br_index::sparse_bitvector_t(first_pos.cbegin(),first_pos.cend());
        }
        fclose(file_ssa_rev);

        // read .rev.esa
        FILE* file_esa_rev = open_aux_file(input_rev.c_str(),EXTESA,"rb");

        idx.samples_lastR = sdsl::int_vector<>(rR,0,log_n);
        idx.last_to_runR = sdsl::int_vector<>(rR,0,log_rR);
        ulint last_SA_valR; // temporary
        {
            auto pos_run_pairs = std::vector<std::pair<ulint,ulint>>(rR);
            /*ulint* sa_val;
            for (ulint i = 0; i < rR; ++i) {
                fread(sa_val,SABYTES,1,file_esa_rev);
                size_t s = fread(sa_val,SABYTES,1,file_esa_rev);
                if (s!=1) die(".rev.esa read failed");
                idx.samples_lastR[i] = *sa_val > 0 ? *sa_val-1 : size-1;
                pos_run_pairs[i] = {*sa_val > 0 ? *sa_val-1 : size-1,i};
            }
            last_SA_valR = *sa_val;*/
            ulint* buf = new ulint[2];
            for (ulint i = 0; i < rR; ++i) {
                size_t s = fread(buf, SABYTES, 2, file_esa_rev);
                if (s!=2) die(".rev.ssa read failed.");
                ulint sa_val = buf[1] > 0 ? buf[1] - 1 : size-1;
                idx.samples_lastR[i] = sa_val;
                pos_run_pairs[i] = {sa_val,i};
            }
            last_SA_valR = buf[1];
            delete[] buf;

            std::sort(pos_run_pairs.begin(),pos_run_pairs.end());

            std::vector<ulint> last_pos;
            for (ulint i = 0; i < rR; ++i) {
                last_pos.push_back(pos_run_pairs[i].first);
                idx.last_to_runR[i] = pos_run_pairs[i].second;
            }
            idx.lastR = br_index::sparse_bitvector_t(last_pos.cbegin(),last_pos.cend());
        }
        fclose(file_esa_rev);

        std::cout << "done.\nBuilding run-length compressed PLCP^R ... " << std::flush;

        // construct run-length encoded PLCP^R
        {
            std::vector<ulint> ones, zeros;

            ulint i = idx.FLR(0); // ISA^R[0]
            assert(i>0);
            uchar c = idx.F_at(i); // T[0]
            ulint i0 = i-1;
            uchar c0 = idx.F_at(i-1); // T[\phi(0)]
            ulint l = 0, prev_l = 0;
            ulint acc0 = 0, acc1 = 0;
            ones.push_back(0);
            for (ulint j = 0; j < size-1; ++j) {
                while (c == c0) {
                    // compute FLR(i) & FLR(i0), compare BWT^R chars
                    l++;
                    ulint p = i - idx.F[c];
                    i = idx.bwtR.select(p,c);
                    c = idx.F_at(i);
                    ulint p0 = i0 - idx.F[c0];
                    i0 = idx.bwtR.select(p0,c0);
                    c0 = idx.F_at(i0);
                }
                if (l==0) {
                    ulint p = i - idx.F[c];
                    i = idx.bwtR.select(p,c); // i = FL(i)
                    c = idx.F_at(i);
                    assert(i>0 || j == size-2);
                    i0 = i - 1;
                    c0 = idx.F_at(i0);
                }

                if (j == 0) {
                    if (l) {
                        acc0 += l - 1;
                        zeros.push_back(acc0);
                        acc1++;
                    } else {
                        acc1++;
                    }
                }
                else if (l + 1 - prev_l) {
                    ones.push_back(acc1);
                    acc0 += l + 1 - prev_l;
                    zeros.push_back(acc0);
                    acc1++;
                } else {
                    acc1++;
                }

                prev_l = l;
                if (l) l--;      
            }
            ones.push_back(acc1);

            idx.plcpR = permuted_lcp<>(size,ones,zeros);
        }

        std::cout << "done.\nBuilding kmer^R[0,bl) ... " << std::flush;


        // consruct kmerR[0,bl) (contraction shortcut bitvector)
        {
            auto kmer_pos = std::vector<std::vector<ulint>>(bl);
            // first run head
            for (ulint i = 0; i < bl; ++i) kmer_pos[i].push_back(0);

            ulint prev_run_start = 0;
            ulint run_start, run_end;
            ulint cur_sa, min_lcp_l, min_lcp_r;

            for (ulint i = 1; i < rR; ++i) {
                auto run_range = idx.bwtR.run_range(i);
                run_start = run_range.first;
                run_end = run_range.second;

                ulint start_sa = idx.samples_firstR[i] < size-1 ? idx.samples_firstR[i]+1 : 0;
                min_lcp_l = min_lcp_r = idx.plcpR[start_sa];
                if (min_lcp_l > bl) min_lcp_l = min_lcp_r = bl;
                
                for (ulint l = min_lcp_l; l < bl; ++l) kmer_pos[l].push_back(run_start);

                // backward scan from run head
                cur_sa = start_sa;
                for (ulint p = run_start-1; p > prev_run_start; --p) {
                    if (min_lcp_l==0) break;

                    cur_sa = idx.PhiR(cur_sa);
                    ulint lcp = idx.plcpR[cur_sa];
                    if (lcp < min_lcp_l) {
                        for (ulint l = lcp; l < min_lcp_l; ++l) {
                            if (kmer_pos[l].back() != p) kmer_pos[l].push_back(p);
                        } 
                        min_lcp_l = lcp;
                    }
                }
                // forward scan from run head
                cur_sa = start_sa;
                for (ulint p = run_start+1; p <= run_end; ++p) {
                    if (min_lcp_r==0) break;

                    cur_sa = idx.PhiIR(cur_sa);
                    ulint lcp = idx.plcpR[cur_sa];
                    if (lcp < min_lcp_r) {
                        for (ulint l = lcp; l < min_lcp_r; ++l) kmer_pos[l].push_back(p);
                        min_lcp_r = lcp;
                    }
                }

                prev_run_start = run_start;
            }
            // last run tail
            for (ulint i = 0; i < bl; ++i) kmer_pos[i].push_back(size);

            idx.kmerR = std::vector<br_index::sparse_bitvector_t>(bl);
            for (ulint i = 0; i < bl; ++i) {
                idx.kmerR[i] = br_index::sparse_bitvector_t(kmer_pos[i].cbegin(),kmer_pos[i].cend());
            }
        }

        std::cout << "done." << std::endl;
        std::cout << "Completed br-index construction." << std::endl;
        
        return *this;
    }

    ulint save_to_file(std::string const& output) {
        std::cout << "Saving br-index to " << output + EXTIDX << " ... " << std::flush;
        std::ofstream f(output + EXTIDX);
        ulint bytes = idx.serialize(f);
        std::cout << "done.\nTotal index size: " << bytes << " bytes." << std::endl;
        return bytes;
    }
};

}; // namespace bri

void print_help(char** argv, Args &args) {
    std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
    std::cout << "Build br-index from Prefix-Free Parsed files." << std::endl;
    std::cout << " .bwt, .ssa, .esa, .rev.bwt, .rev.ssa, .rev.esa are necessary" << std::endl;
    std::cout << " unless you use -i option (in-memory construction, slow&memory consuming)." << std::endl << std::endl;
    std::cout << "  Options: " << std::endl
        << "\t-h  \tshow help and exit" << std::endl
        << "\t-i  \tin-memory construction using constructor of br_index" << std::endl
        << "\t-l L\tparameter bl for contraction shortcut, def. " << args.bl << std::endl
        << "\t-o O\tspecified output file basename, def. <input filename> " << std::endl;
    exit(1);
}

void parse_args( int argc, char** argv, Args& arg ) {
    int c;
    extern char *optarg;
    extern int optind;

    bool specify_out = false;

    puts("==== Command line:");
    for(int i=0;i<argc;i++)
        printf(" %s",argv[i]);
    puts("");

    std::string sarg;
    while ((c = getopt( argc, argv, "hil:o:v") ) != -1) {
        switch(c) {
            case 'l':
            sarg.assign( optarg );
            arg.bl = stoi( sarg ); break;
            case 'o':
            specify_out = true;
            arg.output_base.assign( optarg ); break;
            case 'v':
            arg.verbose++; break;
            case 'h':
            print_help(argv, arg); exit(1);
            case 'i':
            arg.inmemory = true; break;
            case '?':
            std::cout << "Unknown option. Use -h for help." << std::endl;
            exit(1);
        }
    }
    // the only input parameter is the file name 
    if (argc == optind+1) {
        arg.input_file.assign( argv[optind] );
        if (!specify_out) arg.output_base.assign( argv[optind] );
    }
    else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_help(argv,arg);
    }
    // check algorithm parameters 
    if(arg.bl < 0) {
        std::cout << "bl must be nonnegative integer\n";
        exit(1);
    }
}

int main(int argc, char** argv) {
    Args arg;
    parse_args(argc, argv, arg);

    if (arg.inmemory) {
        std::cout << "In-memory construction of br_index ..." << std::endl;

        std::string input;
        std::ifstream fs(arg.input_file);
        std::stringstream buffer;
        buffer << fs.rdbuf();
        input = buffer.str();

        br_index idx(input,arg.bl,false);

        std::cout << "Saving br-index to " << arg.output_base + EXTIDX << " ... " << std::flush;
        std::ofstream f(arg.output_base + "." + EXTIDX);
        ulint bytes = idx.serialize(f);
        std::cout << "done.\nTotal index size: " << bytes << " bytes." << std::endl;

        return 0;
    }

    br_index_builder builder;
    ulint idx_size = builder.build_from_pfp(arg.input_file,arg.bl).save_to_file(arg.output_base);
}
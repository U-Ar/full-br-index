#include "definitions.hpp"
#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include "permuted_lcp.hpp"
#include "utils.hpp"
#include "br_index.hpp"
extern "C" {
    #include "pfpbuilder/utils.h"
}

using namespace bri;

// index file extension
const std::string EXTIDX = "bri";

// struct containing command line parameters and other globals
struct Args {
    std::string input_file = "";
    std::string output_base = "";
    int bl = 8;
    int verbose=0;
    FILE *tmp_parse_file, *last_file, *sa_file; 
};

class br_index_builder {
    br_index idx;
public:
    // input: file basename of Prefix-Free Parsing (.bwt, .ssa, .esa are necessary)
    // bl: parameter
    br_index_builder& build_from_pfp(std::string const& input, int bl) {
        idx = br_index();


        // read .bwt
        std::ifstream fbwt(input + ".bwt", std::ifstream::ate);
        if (!fbwt.is_open()) {
            perror(__func__);
            throw new std::runtime_error("Cannot open file " + input + ".bwt");
        }
        long size = fbwt.tellg();
        fbwt.seekg(0);

        // construct alphabet remapper (null character \0 is mapped to \1)
        idx.remap = std::vector<uchar>(256,0); idx.remap[0] = 1;
        idx.remap_inv = std::vector<uchar>(256,0); idx.remap_inv[1] = 0;
        idx.sigma = 1;
        {
            std::vector<ulint> freqs(256,0); freqs[0] = 1;
            for (size_t i = 0; i < size; ++i) {
                if (freqs[(uchar)fbwt.get()]++ == 0) idx.sigma++;
                if (idx.sigma >= 255) {
                    die("Alphabet overflow");
                }
            }
            uchar new_c = 2;
            for (ulint c = 1; c < 256; ++c) {
                if (freqs[(uchar)c] != 0) {
                    idx.remap[(uchar)c] = new_c;
                    idx.remap_inv[new_c++] = (uchar)c;
                }
            }
        }

        // build RLBWT with remapper
        fbwt.seekg(0);
        idx.bwt.load_from_plain(fbwt,size,idx.remap);

        // read .ssa & .esa
        FILE* file_ssa = open_aux_file(arg.input_file,EXTSSA,"rb");
        FILE* file_esa = open_aux_file(arg.input_file,EXTESA,"rb");



        // load from reversed files
        std::string input_rev = input + ".rev";





        
        return this;
    }

    ulint save_to_file(std::string const& output) {
        std::ofstream f(output + EXTIDX);
        return idx.serialize(f);
    }
};

void print_help(char** argv, Args &args) {
    std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
    std::cout << "  Options: " << std::endl
        << "\t-h  \tshow help and exit" << std::endl
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
    while ((c = getopt( argc, argv, "hl:o:v") ) != -1) {
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
            case '?':
            cout << "Unknown option. Use -h for help." << endl;
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
        std::cout << "bl must be nonnegative value\n";
        exit(1);
    }
}

int main(int argc, char** argv) {
    Args arg;
    parse_args(argc, argv, arg);
    std::cout << "bl: " << arg.bl << std::endl;
    std::cout << "Index file will be saved to " << arg.output_base << "." << EXTIDX << std::endl;

    br_index_builder builder;
    ulint idx_size = builder.build_from_pfp(arg.input_file,arg.bl).save_to_file(arg.output_base);
}
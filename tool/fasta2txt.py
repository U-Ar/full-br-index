import sys



if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("---- FASTA Concatenator ----")
        print("concatenate documents in .fasta file")

        print()
        print("Usage: python fasta2txt.py <input> <output> [<delim>]")
        print("    <input>  input fasta file")
        print("    <output> output file name")
        print("    <delim>  delimiter (none by default)")

        sys.exit(0)
    
    print()

    ifname = sys.argv[1]
    ofname = sys.argv[2]
    delim = ""
    if len(sys.argv) > 3:
        delim = sys.argv[3]

    items = 0

    with open(ifname, mode="r") as fi:
        print("reading {}".format(ifname))
        with open(ofname, mode="w") as fo:
            print("opening {}".format(ofname))
            line = fi.readline()

            while line != "":
                if line[0] == ">":
                    items += 1
                    if items != 1:
                        fo.write(delim)
                else:
                    l = line.strip()
                    if len(l) > 0:
                        fo.write(l)
                line = fi.readline()
    
    print()
    print("successfully processed")
    print("concatenated items: {}".format(items))
    print()
    
    print("----------------------------")
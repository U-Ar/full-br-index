import sys

def print_help():
    print("-- Verifier for bri-count&bri-locate --")
    print("Counts the number of occs with 0, 1, 2 mismatches")
    print("Note: VERY SLOW, use only for verification")
    print()
    print("Usage: python linear_scan.py <text> <patt>")
    print("    <text>  text file")
    print("    <patt>  pattern file in Pizza&Chili format")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    textfile = sys.argv[1]
    pattfile = sys.argv[2]

    print("read {} and {}".format(textfile,pattfile))

    text = ""
    with open(textfile,"r") as f:
        text = f.read()
    
    occ0 = 0; occ1 = 0; occ2 = 0

    with open(pattfile,"r") as f:
        header = f.readline()
        start_pos = header.find("number=")+7
        end_pos = header[start_pos:].find(" ") + start_pos
        n = int(header[start_pos:end_pos])

        start_pos = header.find("length=")+7
        end_pos = header[start_pos:].find(" ") + start_pos
        m = int(header[start_pos:end_pos])

        for i in range(n):
            p = f.read(m)
            o0 = 0; o1 = 0; o2 = 0

            for j in range(len(text)-m+1):
                errors = 0
                for k in range(m):
                    if p[k] != text[j+k]:
                        errors += 1
                if errors == 0:
                    o0 += 1
                elif errors == 1:
                    o1 += 1
                elif errors == 2:
                    o2 += 1
                    
            print("{}-th pattern:".format(i+1))
            print(" - occs with 0 miss : {}".format(o0))
            print(" - occs with 1 miss : {}".format(o1))
            print(" - occs with 2 miss : {}".format(o2))

            occ0 += o0 
            occ1 += o1 
            occ2 += o2


    print("<occs with 0 miss>: {}".format(occ0))
    print("<occs with 1 miss>: {}".format(occ1))
    print("<occs with 2 miss>: {}".format(occ2))
    print("<total occs>      : {}".format(occ0+occ1+occ2))
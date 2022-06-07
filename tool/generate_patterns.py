import sys
import random

def print_help():
    print("-- Pizza&Chili format pattern file generator --")
    print("Randomly extracts substrings from a text file, and save <text>.patt")
    print("Encoding is ISO-8859-1 by default. Edit it as you need.")
    print()
    print("Usage: python generate_patterns.py <text> <n> <m> [<patt>]")
    print("    <text>  text file")
    print("    <n>     number of generated patterns")
    print("    <m>     length of patterns")
    print("    <patt>  output filename (<text>.patt by default)")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    textfile = sys.argv[1]
    n = int(sys.argv[2])
    m = int(sys.argv[3])
    
    outputfile = textfile + ".patt"
    if len(sys.argv) == 5:
        outputfile = sys.argv[4]

    print("read {}".format(textfile))

    text = ""
    with open(textfile,"r", encoding = "ISO-8859-1") as f:
        text = f.read()

    print("write to {}".format(outputfile))

    with open(outputfile,"w", encoding = "ISO-8859-1") as f:
        f.write("# number={} length={} file={}\n".format(n,m,textfile))
        for i in range(n):
            r = random.randint(0,len(text)-m)
            f.write(text[r:r+m])
    
    print("done.")
    
    

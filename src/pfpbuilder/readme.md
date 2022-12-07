# Big-BWT for br-index

Big-BWT(https://gitlab.com/manzai/Big-BWT) modified for the construction of the br-index.

Big-BWT is the tool to build the BWT and optionally the Suffix Array for highly repetitive files using the approach described in *Prefix-Free Parsing for Building Big BWTs* by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini [1]. If you need the detail of PFP refer to the link above.

## Additional features to the original

- bigbwt produces .rev.bwt & .rev.ssa & .rev.esa instead of .bwt & .ssa & .esa if -r option is used. .rev.* files correspond to .* files on the reversed input T^R.
	- pscan.cpp and newscan.cpp reads the input file backwards if -r option is given.
	- BWT checking using simplebwt.c is also implemented for the reversed direction.

## References

\[1\]  Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini 
*Prefix-Free Parsing for Building Big BWTs* [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/).

\[2\] Nong, G., 
*Practical linear-time O(1)-workspace suffix sorting for constant alphabets*, ACM Trans. Inform. Syst., vol. 31, no. 3, pp. 1-15, 2013

HBK_permute is a program that is part of the analysis package of libsequence:

http://molpopgen.org/software/analysis/

This is a minor modification of the code that allows one to use Hudson2001 "polytables" as an input format. An example:

6 6
679 1004 1153 1155 1277 1295 
N N N N N N 
Sim2 A C C G A A
Sim3 G C T G A C
Sim4 A C C T G C
Sim5 A C C T G C
Sim7 A T C T A C
Sim8 A C C G A A

From the libsequence documentationi (http://molpopgen.org/software/libsequence/doc/html/index.html):

"The two numbers on the first line are the sample size and number of segregating sites, respectively. The next line contains the positions of each site. The third line contains the outgroup state for each variable site--if this is unknown, use an 'N' to indicate the ambiguity. The rest of the lines contain a unique sequence name, and then the states of each segregating site."

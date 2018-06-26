import os,sys
f = open(sys.argv[1], "r")
#f = open("../anno/get_ltx_test/test.bed", "r")
CHR = []
START = []
END = []
NAME = []
STRAND = []
LEN = []
ltx = []
for lines in f:
    chr, start, end, name, non, strand = lines.strip().split("\t")
    CHR.append(chr)
    START.append(start)
    END.append(end)
    NAME.append(name)
    length = int(end) - int(start)
    STRAND.append(strand)
    START = map(int, START)
    END = map(int, END)
    LEN.append(length)
f.close()
max = 0
# print(len(NAME)) 10
for i in range(len(NAME)):
    if i == 0:
        max = LEN[i]
    if (i < len(NAME)-1) and (NAME[i] == NAME[i+1]):
        if LEN[i] >= max:
            max = LEN[i]
            indx = i
        else:
            max = max
            indx = indx
    else:
        if LEN[i] >= max:
            max = LEN[i]
            indx = i
        ltx.append(indx)
        max = 0
        index = 0
        continue

oput = open(sys.argv[2],"w")
#oput = open("test_out","w")
for idx in ltx:
    LTX = "\t".join([CHR[idx], str(START[idx]), str(END[idx]), NAME[idx], "0", STRAND[idx]]) +"\n"
    oput.write(LTX)
oput.close()





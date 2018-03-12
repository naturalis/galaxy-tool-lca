#!/usr/bin/python
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input_file', metavar='galaxy blast output', dest='input', type=str,
			help='input data in galaxy blast format', default='', required=True)
parser.add_argument('-o', '--output_file', metavar='output file', dest='output', type=str,
			help='results file in tabular', required=True)
parser.add_argument('-b', '--bitscore', metavar='bitscore top percentage treshold', dest='top', type=str,
			help='top hits to find the lowest common ancestor', required=True)
parser.add_argument('-id', metavar='identity', dest='id', type=str,
			help='identity treshold', required=True)
args = parser.parse_args()

def get_lca(otu):
    #find highest bitscore and calculate lowest bitscore treshold
    highestScore = 0
    for x in otu:
        bitscore = float(x[7])
        highestScore = bitscore if bitscore > highestScore else highestScore
    topTreshold = float(highestScore) * (1 - (float(args.top)/100))

    #place the taxon comlumn in lists for the zip funtion
    taxons = []
    for tax in otu:
        if float(tax[7]) >= topTreshold and float(tax[3]) >= float(args.id):
            taxons.append(map(str.strip, tax[-1].split("/")))
    #use zip function for all* taxon lists
    zipper = zip(*taxons)

    #find the uniq taxon over all levels, if there are more then one uniq taxon a higher taxon is checked
    count = 0
    for y in zipper:
        if len(set(y)) > 1:
            break
        count += 1
    taxonomy = map(str.strip, otu[0][-1].split("/"))[:count]
    taxonLevels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    with open(args.output, "a") as output:
        if taxonomy:
            output.write(otu[0][0]+"\t"+taxonLevels[len(taxonomy)-1]+"\t"+taxonomy[-1]+"\t"+"/".join(taxonomy)+"\n")
        else:
            output.write(otu[0][0] + "\t" + "no identification" + "\t" + "" + "\t" + "" + "\n")

def linecount():
    i = 0
    with open(args.input, "r") as f:
        for i, l in enumerate(f):
            pass
    return i

def lca():
    lastLineCount = linecount()
    with open(args.input, "r") as input:
        otuList = []
        otuLines = []
        for num, line in enumerate(input):
            if line.split("\t")[0] not in otuList and num != 0 or num == lastLineCount:

                if num == lastLineCount:
                    otuList.append(line.split("\t")[0])
                    otuLines.append(line.split("\t"))
                #do stuff with the otu 'block'
                get_lca(otuLines)
                otuList = []
                otuLines = []
                otuList.append(line.split("\t")[0])
                otuLines.append(line.split("\t"))
            else:
                otuList.append(line.split("\t")[0])
                otuLines.append(line.split("\t"))





def main():
    lca()

if __name__ == "__main__":
    main()

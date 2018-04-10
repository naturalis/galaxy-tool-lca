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
parser.add_argument('-t','--tophit', metavar='tophit', dest='tophit', type=str,
			help='Check de best hit first, if it is above the gives treshold the tophit will become the output', required=False, choices=['no', 'yes'], nargs='?', default='no')
parser.add_argument('-tid', metavar='top_hit_identity', dest='topid', type=str,
			help='identity treshold for the tophit', required=False, default='100')
parser.add_argument('-tcov', metavar='top_hit_coverage', dest='topcoverage', type=str,
			help='query coverage treshold for the tophit', required=False,  default='100')

args = parser.parse_args()

def get_lca(otu):
    #find highest bitscore and calculate lowest bitscore treshold
    highestScore = 0
    for x in otu:
        bitscore = float(x[7])
        highestScore = bitscore if bitscore > highestScore else highestScore
    topTreshold = float(highestScore) * (1 - (float(args.top)/100))

    #place the taxon column in lists for the zip function
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
    if taxonomy:
        return otu[0][0] + "\t" + taxonLevels[len(taxonomy) - 1] + "\t" + taxonomy[-1] + "\t" + "/".join(taxonomy) + "\tlca\n"
    else:
        return otu[0][0] + "\t" + "no identification" + "\t" + "" + "\t" + "" + "\tno identification\n"

def check_best_hit(otu):
    if float(otu[0][3]) >= float(args.topid) and float(otu[0][5]) >= float(args.topcoverage):
        taxonomy = map(str.strip, otu[0][-1].split("/"))
        return otu[0][0] + "\tspecies\t" + taxonomy[-1] + "\t" + "/".join(taxonomy) + "\tbest hit\n"
    else:
        return False

def determine_taxonomy(otu):
    with open(args.output, "a") as output:
        if args.tophit == "yes":
            bestHit = check_best_hit(otu)
            if bestHit:
                output.write(bestHit)
            else:
                resultingTaxonomy = get_lca(otu)
                output.write(resultingTaxonomy)
        else:
            resultingTaxonomy = get_lca(otu)
            output.write(resultingTaxonomy)

def linecount():
    i = 0
    with open(args.input, "r") as f:
        for i, l in enumerate(f):
            pass
    return i

def lca():
    """
    This method loops trough the BLAST output and all the hits per otu will be the the input for the get_lca method.
    The first line starts with "Query ID", this is the header so it will not be used. Every line is stored in the
    otuLines variable.
    """
    lastLineCount = linecount()
    with open(args.input, "r") as input:
        otuList = []
        otuLines = []
        for num, line in enumerate(input):
            if line.split("\t")[0].strip() != "Query ID":
                if line.split("\t")[0] not in otuList and num != 1 or num == lastLineCount:
                    if num == lastLineCount:
                        otuList.append(line.split("\t")[0])
                        otuLines.append(line.split("\t"))
                    #do stuff with the otu 'block'
                    determine_taxonomy(otuLines)
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

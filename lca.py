#!/usr/bin/env python

# Define a type to define limits for lca_id_top parameter
def lca_id_top_type(topn: str or int):
    topn = int(topn)
    if topn < 1:
        raise argparse.ArgumentTypeError("Minimum top hits is 1")
    return topn

# Define a type to define limits for lca_id_delta parameter
def lca_id_delta_type(delta: str or float):
    delta = float(delta)
    if delta < 0.0:
        raise argparse.ArgumentTypeError("Minimum top delta is 0.0")
    return delta

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
parser.add_argument('-cov', metavar='coverage', dest='cov', type=str,
            help='coverage treshold', required=True)
parser.add_argument('-t','--tophit', metavar='tophit', dest='tophit', type=str,
            help='Check de best hit first, if it is above the gives treshold the tophit will become the output', required=False, choices=['only_lca', 'lca_threshold', 'best_hit', "best_hits_range"], nargs='?', default='only_lca')
parser.add_argument('-tid', metavar='top_hit_identity', dest='topid', type=str,
            help='identity treshold for the tophit', required=False, default='100')
parser.add_argument('--lca_id_top', type = lca_id_top_type, default = 1, required = False,
            help = 'When using `--tophit lca_threshold`, use this many unique top hit values')
parser.add_argument('--lca_id_delta', type = lca_id_delta_type, default = 0, required = False,
            help = 'When using `--tophit lca_threshold`, allow this amount of deviation below the top hit')
parser.add_argument('-tcov', metavar='top_hit_coverage', dest='topcoverage', type=str,
            help='query coverage treshold for the tophit', required=False,  default='100')
parser.add_argument('-fh', metavar='filter hits', dest='filterHitsParam', type=str,
            help='filter out hit that contain unwanted taxonomy', required=False, default="",nargs='?')
parser.add_argument('-flh', metavar='filter lca hits', dest='filterLcaHits', type=str,
            help='do not use a String in de lca determination', required=False, default="",nargs='?')
parser.add_argument('-fs', metavar='filter on taxonomy source', dest='filterSourceHits', type=str,
            help='do not use hit when taxonomy from source', required=False, default="",nargs='?')
parser.add_argument('-minbit', dest='minbit', type=str, required=False, nargs='?', default="0")

args = parser.parse_args()

def filter_check(filterParam, line):
    """
    This method checks if the input (line) does not contain a certain string that is present
    in filterParam. filterParam is a comma separated string, it will be converted to a list by
    splitting it on the "," character. If a string present in filterParam is also present in line
    it will return False.
    """
    filterHitsParam = filterParam.strip()
    if filterParam[-1] == ",":
        filterHitsParam = filterParam[:-1]
    filter = filterHitsParam.split(",")
    for x in filter:
        if x.lower() in line.lower():
            a = False
            break
        else:
            a = True
    return a


def remove_hits(otu):
    """
    This method loops through the lines of one otu. The taxonomy column will be the input for the filter_check() method.
    If the method filter_check() returns True the line will be appended to the filteredOtu list.
    """
    filteredOtu = []
    for line in otu:
        if filter_check(args.filterHitsParam.strip(), line[-1]):
            filteredOtu.append(line)
    return filteredOtu


def remove_wrong_source_hits(otu):
    filteredOtu = []
    for line in otu:
        if filter_check(args.filterSourceHits, line[8]):
            filteredOtu.append(line)
    return filteredOtu


def remove_taxon(zippedTaxonomy):
    """

    """
    filteredZipper = []
    for level in zippedTaxonomy:
        filteredRank = []
        for item in level:
            if filter_check(args.filterLcaHits, item):
                filteredRank.append(item)
        if not filteredRank:
            filteredRank.append("no identification")
        filteredZipper.append(filteredRank)
    return filteredZipper


def check_best_hit(otu):
    """
    This method will be used if the user chooses parameter -t best_hit.
    Before checking the LCA there will be checked if the top hit passes
    some thresholds (args.topid, args.topcoverage). If yes, the top hit will be returned
    and gets the tag "best hit" in the last column.
    """
    if float(otu[0][4]) >= float(args.topid) and float(otu[0][5]) >= float(args.topcoverage):
        taxonomy = list(map(str.strip, otu[0][-1].split(" / ")))
        return otu[0][0] + "\tspecies\t" + taxonomy[-1] + "\t" + "\t".join(taxonomy) + "\tbest hit\n"
    else:
        return False


def check_best_hit_range(otu_filtered):
    """
    This method will be used if the user chooses -t best_hit_range.
    All the hits will be checked if they pass a certain threshold.
    All the hits the pass the threshold with the same taxonomy will be merged to one line.
    The lowest and highest scores (range) will be written to the last column.
    """
    hitList = {}
    output = ""
    for x in otu_filtered:
        if float(x[4]) >= float(args.topid) and float(x[5]) >= float(args.topcoverage):
            taxonomy = list(map(str.strip, x[-1].split(" / ")))
            # for x in otu_filtered:#loop trough otus
            #     if float(x[4]) >= float(args.topid) and float(x[5]) >= float(args.topcoverage):#check thresholds
            #         taxonomy = map(str.strip, x[-1].split(" / "))
            line = x[0] + "\tspecies\t" + taxonomy[-1] + "\t" + "\t".join(taxonomy) + "\ttop hit"
            if line not in hitList:  # check if the taxonomy was already present before, if no add, if yes append
                hitList[line] = [[float(x[4])], [float(x[5])]]
            else:
                hitList[line][0].append(float(x[4]))
                hitList[line][1].append(float(x[5]))
    for y in hitList:  # loop trough otus again, this time a few will be merged by taxonomy
        numberOfHits = str(len(hitList[y][0]))
        best = y.replace("top hit", "top hit (" + numberOfHits + ")")
        # put the lowest and highest scores to the output
        line = best + "\t" + str(round(min(hitList[y][0]), 1)) + "-" + str(round(max(hitList[y][0]), 1)) + "\t" + str(
            round(min(hitList[y][1]), 1)) + "-" + str(round(max(hitList[y][1]), 1)) + "\n"
        output += line
    if output:
        return output
    else:
        return False


def calculate_bitscore_threshold(otu):
    highestScore = 0
    for x in otu:
        bitscore = float(x[7])
        highestScore = bitscore if bitscore > highestScore else highestScore
    topthreshold = float(highestScore) * (1 - (float(args.top) / 100))
    return topthreshold


def zip_taxonomy_column(otu, topthreshold):
    """
    Of all the otus that passes the threshold the taxonomy will be "zipped" https://www.w3schools.com/python/ref_func_zip.asp
    A zipped list (zippedTaxonomy) will look like something like this: [[family, family], [genus, genus], [species, species]]
    """
    taxons = []
    for tax in otu:
        # if float(tax[7]) >= topthreshold and float(tax[4]) >= float(args.id) and float(tax[5]) >= float(args.cov) and float(tax[7]) >= 50:#float(args.minbit):
        if float(tax[7]) >= topthreshold and float(tax[4]) >= float(args.id) and float(tax[5]) >= float(
                args.cov) and float(tax[7]) >= float(args.minbit):
            taxons.append(list(map(str.strip, tax[-1].split(" / "))))
        # if float(tax[7]) >= topthreshold and float(tax[4]) >= float(args.id) and float(tax[5]) >= float(args.cov) and float(tax[7]) >= float(args.minbit):
        #     taxons.append(map(str.strip, tax[-1].split(" / ")))
    # use zip function for all* taxon lists
    zippedTaxonomy = list(zip(*taxons))
    return zippedTaxonomy


def find_lca(zippedTaxonomy):
    """
    This method finds the lowest common ancestor. It make the lists in zippedTaxonomy unique
    after that if there are more then one different taxonomies left the loop stops.
    example with a zipped file like this:[[family1, family1], [genus1, genus2], [species1, species2]]
    When you make family unique the will be only one family left (family1), if you make genus unique
    there will be two genera left. So the lowest common ancesor is family1.
    """
    count = 0
    taxonomy = []
    if zippedTaxonomy:
        # change species to no identification, lca starts from genus
        zippedTaxonomy[-1] = "no identification"
    for y in zippedTaxonomy:
        # make unique and check if there are more then 1
        if len(set(y)) > 1:  # if more then 1 stop checking
            break
        count += 1
        taxonomy.append(list(set(y))[0])  # put the lca taxonomy in a variable
    t = True
    newTaxonomy = []
    for rev in reversed(taxonomy):
        if rev == "no identification" and t:
            pass
        else:
            t = False
            newTaxonomy.append(rev)
    return list(reversed(newTaxonomy))


def generate_output_line(taxonomy, otu):
    taxonLevels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxonLevel = len(taxonomy) - 1
    while len(taxonomy) < 7:
        taxonomy.append("no identification")

    if all(rank == "no identification" for rank in taxonomy):
        return otu[0][0] + "\tno identification\tno identification\t" + "\t".join(taxonomy).strip() + "\tno lca"
    else:
        if taxonLevel == 6:
            taxonomy[-1] = ""
            return otu[0][0] + "\t" + taxonLevels[taxonLevel - 1] + "\t" + taxonomy[taxonLevel - 1] + "\t" + "\t".join(
                taxonomy).strip() + "\tno identification\tlca"
        else:
            return otu[0][0] + "\t" + taxonLevels[taxonLevel] + "\t" + taxonomy[taxonLevel] + "\t" + "\t".join(
                taxonomy) + "\tlca"


def get_lca(otu):
    """
    This method contains a few other methods to determine the lca.
    """
    # find highest bitscore and calculate lowest bitscore threshold
    topthreshold = calculate_bitscore_threshold(otu)

    # place the taxon column in lists for the zip function
    zippedTaxonomy = zip_taxonomy_column(otu, topthreshold)

    # filter the taxonomy levels, taxons that match 'filter lca hits' parameters are removed
    if args.filterLcaHits:
        zippedTaxonomy = remove_taxon(zippedTaxonomy)
    # find the lca and make the line for the output file
    outputLine = generate_output_line(find_lca(zippedTaxonomy), otu)
    return outputLine

def threshold_filter(otu_array: list[list[str]], identity_array: dict[str, list[float]], topn: int = 1, topdelta: float = 0) -> list[list[str]]:
    """
    Filter OTU hits by the percentage identity found among the hits for individual OTUs. By default, the highest unique percentage identity is used.
    Using the topn argument, users can select a lower value among the given unique percentage identities. The topdelta parameter can be used to subtract
    a set value from the resulting percentage identity threshold.
    Example: given the following array of percentage identities [100.0, 99.0, 98.0, 90.0] the following final threshold will be computed:
        topn = 1, topdelta = 0.0: 100.0
        topn = 3, topdelta = 0.0: 98.0
        topn = 3, topdelta = 20.0: 78.0
    OTU hits below the resulting threshold will not be used to infer taxonomy based on LCA. 
    """
    filtered_otu_array = []
    for otu in otu_array:
        # print(float(min(sorted(set(identity_array[otu[0]]), reverse = True)[0:topn])) - topdelta) # Uncomment this line to check percentage identity value
        if float(otu[4]) >= float(min(sorted(set(identity_array[otu[0]]), reverse = True)[0:topn])) - topdelta:
            filtered_otu_array.append(otu)
    return filtered_otu_array

def determine_taxonomy(otu, identity_array = None):
    """
    This method contains other methods to determine the output.
    The first step is to do some filtering, if a line contains a certain word it will be removed.
    If there are no lines left the otu gets no identification.
    :param otu: list with all the lines of one otu
    """
    otu_filtered = otu
    if args.filterHitsParam.strip() and args.filterHitsParam.strip() != "none":
        otu_filtered = remove_hits(otu_filtered)  # filter on args.filterHitsParam
    if args.filterSourceHits.strip():
        otu_filtered = remove_wrong_source_hits(otu_filtered)  # filter on args.filterSourceHits
    endLine = "\n"
    with open(args.output, "a") as output:
        bestHit = ""
        if otu_filtered:
            if args.tophit == "best_hit":
                bestHit = check_best_hit(otu_filtered)
            elif args.tophit == "best_hits_range":
                endLine = "\t\t\n"
                bestHit = check_best_hit_range(otu_filtered)
            elif args.tophit == "lca_threshold":
                threshold_passed_otu = threshold_filter(otu_filtered, identity_array, args.lca_id_top, args.lca_id_delta)
                if threshold_passed_otu:
                    resultingTaxonomy = get_lca(threshold_passed_otu)
                    output.write(resultingTaxonomy+endLine)
            else:
                resultingTaxonomy = get_lca(otu_filtered)
                output.write(resultingTaxonomy + endLine)
            if bestHit:  # if there was a best hit write the best hit to the output
                output.write(bestHit)
            elif args.tophit == "best_hit" or args.tophit == "best_hits_range":  # if there was no best hit, check the lca
                resultingTaxonomy = get_lca(otu_filtered)
                output.write(resultingTaxonomy + endLine)

        else:
            taxonomy = ["no identification", "no identification", "no identification", "no identification",
                        "no identification", "no identification", "no identification"]
            output.write(
                otu[0][0] + "\tno identification\tno identification\t" + "\t".join(
                    taxonomy).strip() + "\tfiltered out" + endLine)


def linecount():
    """
    This method counts the total amount of line in the input file.
    """
    i = 0
    with open(args.input, "r") as f:
        for i, l in enumerate(f):
            pass
    return i


def write_header():
    """
    Write a header line to the output file
    """
    with open(args.output, "a") as output:
        if args.tophit == "best_hits_range":
            output.write(
                "#Query\t#lca rank\t#lca taxon\t#kingdom\t#phylum\t#class\t#order\t#family\t#genus\t#species\t#method\t#identity\t#coverage\n")
        else:
            output.write(
                "#Query\t#lca rank\t#lca taxon\t#kingdom\t#phylum\t#class\t#order\t#family\t#genus\t#species\t#method\n")


def get_highest_identity_per_otu(blast_otu_file_path: str) -> dict[str, list[float]]:
    """
    This function reads the lines of the BLAST OTU table supplied with the -i parameter. The output will be a dictionary where each key is
    an OTU id found in the input file, and the value will be a list containing the identity percentages from all hits of a specific OTU id.
    """
    identity_per_otu = {}
    with open(blast_otu_file_path) as input:
        for line in input:
            if line.split("\t")[0].strip() != "#Query ID":
                line_items = line.split("\t")
                if line_items[0] not in identity_per_otu.keys():
                    identity_per_otu[line_items[0]] = []
                identity_per_otu[line_items[0]].append(float(line_items[4]))
    return identity_per_otu


def get_highest_identity_per_otu(blast_otu_file_path: str) -> dict[str, list[float]]:
    """
    This function reads the lines of the BLAST OTU table supplied with the -i parameter. The output will be a dictionary where each key is
    an OTU id found in the input file, and the value will be a list containing the identity percentages from all hits of a specific OTU id.
    """
    identity_per_otu = {}
    with open(blast_otu_file_path) as input:
        for line in input:
            if line.split("\t")[0].strip() != "#Query ID":
                line_items = line.split("\t")
                if line_items[0] not in identity_per_otu.keys():
                    identity_per_otu[line_items[0]] = []
                identity_per_otu[line_items[0]].append(float(line_items[4]))
    return identity_per_otu


def lca():
    """
    This method loops through the BLAST output and all the hits per otu will be the input for the determine_taxonomy method.
    The first line starts with "Query ID", this is the header so it will not be used. Every line of the same otu is stored in the
    otuLines variable. There are multible otus in one file, only the lines of the same otu need to go in the determine_taxonomy method.
    This works as followed: You loop through the file (blast output) the loop just started so num is not equal to 1 so the else block
    is used and the first line of otu1 and the name 'otu1' will be added to the lists. The second line is still otu1 and otu is present in
    the otuList so again the else block will be used. When the loop reaches otu2, otu 2 is not present in otuLines so the otu1 line go into
    the determine_taxonomy(otuLines) method and the lists will be emptied and the lines of otu2 will now be stored in the list.
    """
    write_header()
    lastLineCount = linecount()
    if args.tophit == "lca_threshold":
        identity_value_list = get_highest_identity_per_otu(args.input)
    with open(args.input, "r") as input:
        otuList = []
        otuLines = []
        for num, line in enumerate(input):  # start loop
            if line.split("\t")[0].strip() != "#Query ID":  # not use the header
                if line.split("\t")[0] not in otuList and num != 1 or num == lastLineCount:
                    if num == lastLineCount:
                        otuList.append(line.split("\t")[0])
                        otuLines.append(line.split("\t"))
                    if args.tophit == "lca_threshold":
                        determine_taxonomy(otuLines, identity_value_list)
                    else:
                        determine_taxonomy(otuLines)#find the lca for the query
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

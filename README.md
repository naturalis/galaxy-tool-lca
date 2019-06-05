# galaxy-tool-lca
A tool to determine the lowest common ancestor from BLAST results with taxonomy. This approach is partly based on MEGAN's (Huson et al., 2007) LCA method. This tool is more flexible and easier to use then MEGAN, it does not specifically use taxonids or separated mapping files. Instead, the script handles files with the taxonomy present in the last column of the input file. Those input files can be generated with our BLAST pipeline (not been published yet) or you can create them manually. Determining the lowest common ancestor can help to identify sequences that do not have a significant good blast hit. The basis of the approach for the analyses goes into the folowing order:
<br /><br />
**Step 1, determine the top percentage based on bitscore:** <br />
Of all the blast hits a sub-selection will be made of the top hits based on a percentage of the bitscore. If the top percentage is set to 8% the top will be calculated like 0.92\*bitscore best hit=top threshold. In this example, all hits with a bitscore above 294.4 will continue to the next step.   

| Query | Subject | Identity percentage | Coverage | bitscore |
| --- | --- | --- | --- | --- |
| Otu1 | hit1 | 93 | 95 | 320 |
| Otu1 | hit2 | 93 | 95 | 320 |
| Otu1 | hit3 | 91 | 94 | 310 |
| Otu1 | hit4 | 90 | 95 | 301 |
| Otu1 | hit5 | 89 | 94 | 300 |
| Otu1 | hit6 | 85 | 95 | 290 |
| Otu1 | hit7 | 84 | 95 | 290 |
| Otu1 | hit8 | 81 | 90 | 281 |
| Otu1 | hit9 | 80 | 89 | 275 |
| Otu1 | hit10 | 79 | 88 | 274 |

<br />

**Step 2, filter on threshold:** <br />
From the top hits there will be filtered on identity, coverage and bitscore. Let's choose a threshold of 90 identity, 90 coverage and 250 bitscore. All hits with scores above the thresholds will go to the next step. 

| Query | Subject | Identity percentage | Coverage | bitscore |
| --- | --- | --- | --- | --- |
| Otu1 | hit1 | 93 | 95 | 320 |
| Otu1 | hit2 | 93 | 95 | 320 |
| Otu1 | hit3 | 91 | 94 | 310 |
| Otu1 | hit4 | 90 | 95 | 301 |
| Otu1 | hit5 | 89 | 94 | 300 |

<br />

**Step 3, determine the lowest common ancestor:** <br />
Of all the remaining hits the lowest common ancestor is determined. The script starts at species level, it checks if all the hits are coming from the same species. If not it checks the genus level and so on. 

| Query | Subject | taxonomy |
| --- | --- | --- |
| Otu1 | hit1 | Animalia / Arthropoda / Insecta / Diptera / Asilidae / Scarbroughia / Scarbroughia delicatula |
| Otu1 | hit2 | Animalia / Arthropoda / Insecta / Diptera / Asilidae / Scarbroughia / Scarbroughia delicatula |
| Otu1 | hit3 | Animalia / Arthropoda / Insecta / Diptera / Asilidae / Schildia / Schildia fragilis |
| Otu1 | hit4 | Animalia / Arthropoda / Insecta / Diptera / Asilidae / Scytomedes / Scytomedes haemorrhoidalis |

With this example data and thresholds the classification of otu1 will be Asilidae. This is a basic explanation of the LCA approach, more advanced filter settings and output options are shown at the example commands.
<br />

## Getting Started
### Prerequisites
python 2.7
### Download
```
git clone https://github.com/naturalis/galaxy-tool-lca
```

## Usage
There is an example input file included in the example folder, this file will be used to execute the example commands. The file consist of 11 columns, has a specific header and is tab separated.<br />
<br />
**Column explanation input file:**

| Column name | Description |
| --- | --- |
| Query ID | Id of the input sequence (qseqid) |
| Subject | means Subject Title (stitle) |
| Subject accession | Subject accession (sacc) |
| Subject Taxonomy ID | unique Subject Taxonomy ID (staxid), this value can be any value and is not used the find the lca |
| Identity percentage | Identity percentage of hit () |
| Coverage | means Query Coverage Per Subject (qcovs) |
| Evalue | expected value of the hit (evalue) |
| Bitscore | Bit score of hit (bitscore) |
| Source | The source of where the taxonomy comes from, this value can be any value and is not used the find the lca |
| Taxonomy | This column contains the taxonomy of the hit in the order Kingdom / phylum / class / order / family / genus /species |

The script itself has multiple parameter options.<br />
<br />
**Parameters:**

| Parameter | Description |
| --- | --- |
| -i | Input file |
| -o | Output file |
| -b | bitscore top percentage threshold |
| -id | Minimum identity threshold |
| -cov | Minimum coverage threshold |
| -t | Check the top hit first or perform an lca analysis on all hits. Options:['only_lca', 'best_hit', "best_hits_range"] |
| -tid | Identity threshold for the tophit, only used when -t best_hit or best_hits_range |
| -tcov | Coverage threshold for the tophit, only used when -t best_hit or best_hits_range |
| -fh | Filter hits, filter out lines that contain a certain string |
| -flh | Filter lca hits, during the determination of the lca ignore this string  |
| -minbit | Minimum bitscore threshold |

### Examples
**Example 1:**<br />
This command performs an lca analysis on all the hits per otu. The top percentage threshold is 8%, minimal idenity is 80% and minimal coverage is set to 80. If you perform an lca analysis on all the hits the classification will never be on species level. 
```
python lca.py -i example/example.tabular -o output1_example.tabular -b 8 -id 80 -cov 80 -t only_lca
```
**Example 1 output explanation:**<br />
Otu 1 has no identification and there is no lca analysis performed. You can see this in the last column. The hits of this otu did not passed the thresholds.

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu1 | no identification | no identification | no identification | no identification | no identification | no identification | no identification | no identification | no identification | no lca |

For Otu3 and Otu 6 an lca analysis is performed. The identification of otu 3 is on genus level and is classified as Chaetopterus. 

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu3 | genus | Chaetopterus | Eukaryota | Annelida | Polychaeta | Spionida | Chaetopteridae | Chaetopterus | no identification | lca |
| Otu6 | class | Polychaeta | Eukaryota | Annelida | Polychaeta | no identification | no identification | no identification | no identification | lca |

Otu16 is classified as "unknown genus". This is because that is how the information is stored in the database. In this case the hit is coming from BOLD http://www.boldsystems.org/index.php/Public_RecordView?processid=POLNB1246-14

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu16 | genus | unknown genus | Animalia | Annelida | Polychaeta | Polychaeta incertae sedis | unknown family | unknown genus | no identification | lca |

**Example 2 taxonomy filter:**<br />
This command performs an lca analysis on all the hits per otu like example 1, but now during the lca analysis if the string "unknown" exist in a rank it will be ignored. In many databases like Genbank or BOLD not all taxonomy is filled in by the uploader. 
```
python lca.py -i example/example.tabular -o output2_example.tabular -b 8 -id 80 -cov 80 -t only_lca -flh unknown 
```
**Example 2 taxonomy filter output explanation:**<br />
Otu6 is now now classified as Thelepus at genus level. If you look at the example.tabular file you can see that this makes sense, in example 1 the rank was pushed up to class level because of the "unknowns" in the first hit. 

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu6 | genus | Thelepus | Eukaryota | Annelida | Polychaeta | Terebellida | Terebellidae | Thelepus | no identification | lca |
| Otu16 | order | Polychaeta incertae sedis | Animalia | Annelida | Polychaeta | Polychaeta incertae sedis | no identification | no identification | no identification | lca |

**Example 3 best hit species level identification:**<br />
This command performs an lca analysis if the top hit falls below the "top hit" thresholds. If the top blast hit (sorted on evalue) exceeds the thresholds the top hit is chosen and the input sequence gets a species level identification. Notice the change in parameters ```-t best_hit -tid 98 -tcov 100```. If the top hit has an identity above 98% en a coverage above or equal to 100% this hit will be the taxonomy classification and will be at species level.    
```
python lca.py -i example/example.tabular -o output3_example.tabular -b 8 -id 80 -cov 80 -t best_hit -tid 98 -tcov 100 -flh unknown 
```
**Example 3 best hit species level identification output explanation:**<br />
Now we have a otu that is identified to species level. Also the last column of Otu9 is different than the previous examples. The value "best hit" means that there is no lca analysis performed, but the species of the top hit is chosen for identification. If you look at the example.tabular file you can see that the top hit of otu9 has an identity above 98 and coverage above or equal to 100. 

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu6 | genus | Thelepus | Eukaryota | Annelida | Polychaeta | Terebellida | Terebellidae | Thelepus | no identification | lca |
| Otu9 | species | Myxine glutinosa | Eukaryota | Chordata | unknown class | Myxiniformes | Myxinidae | Myxine | Myxine glutinosa | best hit |

**Example 4 best hit range:**<br />
In some cases the sequences in the reference database are wrongly morphological identified or contaminated by human or bacterial DNA. It can also happen that a certain marker (16S, CO1, ITS) is not distinctive enough. Many people choose the top hit as the identification for the input sequence without looking at the second hit while then can even be a better choice. In these cases the following command can help to solve this. Notice the changed parameter ```-t best_hits_range```
```
python lca.py -i example/example.tabular -o output4_example.tabular -b 8 -id 80 -cov 80 -t best_hits_range -tid 98 -tcov 100 -flh unknown
```
**Example 4 best hit range output explanation:**<br />
After executing the command otu6 stays the same. Otu6 does not have hits above ```-tid 98 -tcov 100``` so the parameter ```-t best_hits_range``` has no effect here. If we look at otu9 you see that it occurs twice in the output with two extra columns. Otu9 had multiple hits above ```-tid 98 -tcov 100``` on two different species. In the blast output of otu9 there were 6 hits on the species *Myxine glutinosa* with an identity between 98.7 and 100%.   

| #query | #lca rank | #lca taxon | #kingdom | #phylum | #class | #order | #family | #genus | #species | #method | #identity | #coverage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Otu6 | genus | Thelepus | Eukaryota | Annelida | Polychaeta | Terebellida | Terebellidae | Thelepus | no identification | lca |  |  |
| Otu9 | species | Myxine glutinosa | Eukaryota | Chordata | unknown class | Myxiniformes | Myxinidae | Myxine | Myxine glutinosa | top hit (6) | 98.7-100.0 | 100.0-100.0 |
| Otu9 | species | Myxine limosa | Eukaryota | Chordata | unknown class | Myxiniformes | Myxinidae | Myxine | Myxine limosa | top hit (4) | 98.1-98.7 | 100.0-100.0 |

## Author


## How to cite


## Source
Huson, D. H., Auch, A. F., Qi, J., & Schuster, S. C. (2007). MEGAN analysis of metagenomic data. Genome Research, 17(3), 377–386. http://doi.org/10.1101/gr.5969107

## Getting Started
### Prerequisites
This tool uses specifically the BLAST output as input
### Installing
Installing the tool for use in Galaxy
```
cd /home/galaxy/Tools
```
```
git clone https://github.com/naturalis/galaxy-tool-lca
```
Add the following line to /home/galaxy/galaxy/config/tool_conf.xml
```
<tool file="/home/galaxy/Tools/galaxy-tool-lca/lca.xml" />
```

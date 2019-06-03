# galaxy-tool-lca
A tool to determine the lowest common ancestor from BLAST results with taxonomy. This approach is partly based on MEGAN's (Huson et al., 2007) LCA method. This tool is more flexible and easier to use then MEGAN, it does not specificly use taxonids or sepperated mapping files. Instead the script handles files with the taxonomy present in de last column of the input file. Those input files can be generated with our BLAST pipeline (not been published yet) or you can create them manually.   

## Getting Started
### Prerequisites
python 2.7
### Download
git clone https://github.com/naturalis/galaxy-tool-lca

## Usage
### Quick start

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

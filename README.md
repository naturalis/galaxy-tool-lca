# galaxy-tool-lca
A tool to determine the lowest common ancestor from BLAST results with taxonomy. This approach is partly based on MEGAN's (Huson et al., 2007) LCA method. This tool does not specifically use taxonids or separated mapping files. Instead, the script handles files with the taxonomy present in the last column of the input file. Those input files can be generated with our BLAST pipeline or you can create them manually.

## Installation
### Manual
Clone this repo in your Galaxy ***Tools*** directory:  
`git clone https://github.com/naturalis/galaxy-tool-lca`  

Make the scripts executable:  
`chmod 755 galaxy-tool-lca/lca.sh`  
`chmod 755 galaxy-tool-lca/lca.py`  

Append the file ***tool_conf.xml***:    
`<tool file="/path/to/Tools/galaxy-tool-lca/lca.xml" />`  

### Ansible
Depending on your setup the [ansible.builtin.git](https://docs.ansible.com/ansible/latest/collections/ansible/builtin/git_module.html) module could be used.  
[Install the tool](https://docs.ansible.com/ansible/latest/collections/ansible/builtin/git_module.html#examples) by including the following in your dedicated ***.yml** file:  

`  - repo: https://github.com/naturalis/galaxy-tool-lca`  
&ensp;&ensp;`file: lca.xml`  
&ensp;&ensp;`version: master`  

Information on tool usage, parameters and examples can be found [here](https://github.com/naturalis/galaxy-tool-lca/blob/master/usage.md)  

## How to cite
Beentjes KK, Speksnijder AGCL, Schilthuizen M, Hoogeveen M, Pastoor R, et al. (2019) Increased performance of DNA metabarcoding of macroinvertebrates by taxonomic sorting. PLOS ONE 14(12): e0226527. https://doi.org/10.1371/journal.pone.0226527

## Source
Huson, D. H., Auch, A. F., Qi, J., & Schuster, S. C. (2007). MEGAN analysis of metagenomic data. Genome Research, 17(3), 377–386. http://doi.org/10.1101/gr.5969107

# galaxy-tool-lca
A tool to determine the lowest common ancestor from BLAST results
## Getting Started
### Prerequisites
This tool uses specifically the BLAST output as input
### Installing
Installing the tool for use in Galaxy
```
cd /home/galaxy/Tools
```
```
sudo git clone https://github.com/naturalis/galaxy-tool-lca
```
```
sudo chmod 777 galaxy-tool-lca/lca.py
```
```
sudo ln -s /home/galaxy/Tools/galaxy-tool-lca/lca.py /usr/local/bin/lca.py
```
```
sudo cp galaxy-tool-lca/lca.sh /home/galaxy/galaxy/tools/identify/lca.sh
sudo cp galaxy-tool-lca/lca.xml /home/galaxy/galaxy/tools/identify/lca.xml
```
Add the following line to /home/galaxy/galaxy/config/tool_conf.xml
```
<tool file="identify/lca.xml" />
```

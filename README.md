# galaxy-tool-lca
A tool to determine the lowest common ancestor from BLAST results.This repo can be used for the new (03-04-2019) galaxy 19.01 Naturalis server. The old galaxy 16.04 server is not supported anymore with this pipeline.
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
Add the following line to /home/galaxy/galaxy/config/tool_conf.xml
```
<tool file="/home/galaxy/Tools/galaxy-tool-lca/lca.xml" />
```

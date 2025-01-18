# Phylogenetic pipeline

## System requirements

Before starting, make sure the following tools are installed:
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [MCL](https://micans.org/mcl/)

You can install them using:
```bash
# Linux (Debian/Ubuntu)
sudo apt install ncbi-blast+ mcl

# or manually:
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
tar -xzf ncbi-blast-*-x64-linux.tar.gz
export PATH=$PATH:/path/to/blast/bin

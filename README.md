# profun


Test Environment
--------------------------------------------------------------------------------------
Red Hat Enterprise Linux Server release 6.4 (Santiago)

Installation Steps
--------------------------------------------------------------------------------------


**(A) Download and Unzip profun source package**  

Create a working directory called 'profun' where all scripts, programs and databases will reside:

Download the profun code:
```
cd ~/
git clone https://github.com/multicom-toolbox/profun.git
cd profun
```

**(B) Download programs**
```
cd programs
http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/blast%2B/2.2.26/ncbi-blast-2.2.26%2B-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.2.26%2B-x64-linux.tar.gz

```

**(C) Download database**

```
mkdir database
cd database
mkdir swiss_prot
cd swiss_prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz
gzip -d uniprot_sprot.dat.gz
gzip -d uniprot_sprot.fasta.gz

../blast-2.2.26/bin/formatdb  -p T -i uniprot_sprot.fasta

```



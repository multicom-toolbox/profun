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
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
tar -zxvf blast-2.2.26-x64-linux.tar.gz

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

../../programs/blast-2.2.26/bin/formatdb  -p T -i uniprot_sprot.fasta

```

**(D) Test**
```
(1) Run blast
mkdir -p ./test/2SN3-A/sequence ./test/2SN3-A/Blast_output
perl scripts/blast_search_swiss_prot.pl ./test/2SN3-A.fasta  programs/blast-2.2.26/bin/blastpgp ?/test/2SN3-A/sequence /test/2SN3-A/Blast_output ./database/uniprot_sprot/uniprot_sprot.fasta ./database/uniprot_sprot/uniprot_sprot.dat /test/2SN3-A/LOG/BLAST.LOG ./database/uniprot_sprot/uniprot_sprot.fasta   

(2) Run prediction
mkdir ./test/2SN3-A/pred
perl bin/R_all_three_methods.pl ./scripts/ database/ ./test/2SN3-A/sequence ./test/2SN3-A/Blast_output ./test/2SN3-A/pred 

```

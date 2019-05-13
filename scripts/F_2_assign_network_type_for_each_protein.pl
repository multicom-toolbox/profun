#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 use POSIX qw/ceil/; #this is for function ceil()



  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args
      print("This program will load all protein ID, and also load all networks, assign the type of network for each protein ...\n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_uniprot_dat dir_networks(*.network, could be gene-gene interaction network, protein-protein interaction network) protein_list addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../data/data_for_F_2_all_networks ../result/F_1_all_proteins_from_BLAST ../result/F_2_protein_lists_with_network_type\n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";
 
 my($addr_uniprot)=$ARGV[0];
 my($dir_networks)=$ARGV[1];      # all networks
 my($addr_list)=$ARGV[2];         # protein lists for make assignment
 my($addr_output)=$ARGV[3];

 -s $addr_uniprot || die "Cannot open $addr_uniprot!\n";
 -s $dir_networks || die "Cannot open $dir_networks!\n";
 -s $addr_list || die "Cannot open $addr_list!\n";

 
 my($key,$value,$IN,$OUT,$path_input,$target_name,$path_output,$line,$file,$i,$j,$evalue,$network_type);
 my(@tem_split,@tem_2,@files);
 $OUT = new FileHandle ">$addr_output";
 defined($OUT) || die "cannot open output file $addr_output\n";
 $OUT->close(); 

 $| = 1;
 print "Loading uniprot database ... ";
 #sleep(1);
 #############################################################################################
 #
 #                1. Load uniprot database, get hash table to convert AC to ID, ID to geneID, geneID to ID, AC to AC line, geneID to AC line
 #
 #############################################################################################
 my(%AC2ID)=();                # convert protein AC to uniprot ID
 my(%ID2gene)=();              # convert uniprot ID to geneID
 my(%gene2ID)=();              # convert geneID to uniprot ID

 my(%AC2ACline)=();            # convert protein AC to the AC line.
 my(%gene2ACline)=();          # convert geneID to the AC line

 
 my($v_ID,$v_gene,$index_AC,$line_AC);
 my(@v_AC);

 $v_ID="NULL";
 @v_AC=();
 $v_gene="NULL";
 $line_AC="NULL";
 $index_AC=0;

 $IN=new FileHandle "$addr_uniprot";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
     @tem_split=split(/\s+/,$line);
	 if(substr($line,0,2) eq "//")
	 {# new line comes
		 if($v_ID eq "NULL")
		 {
			 print "Severe warning! no ID, check @v_AC,$v_gene \n";
			 exit(0);
		 }
		 if($line_AC eq "NULL")
		 {
			 print "Severe error, no AC ?? @v_AC,$v_gene,$v_ID\n";
			 exit(0);
		 }
         for($i=0;$i<@v_AC;$i++)
		 {# build relationship of AC and ID
			 if(not exists $AC2ID{$v_AC[$i]})
			 {
				 $AC2ID{$v_AC[$i]}=$v_ID;
			 }
			 else
			 {
				 #print "!!!!!! this protein AC $v_AC[$i] has more than one ID !!!!\n";

			 }
			 if(not exists $AC2ACline{$v_AC[$i]})
			 {
				 $AC2ACline{$v_AC[$i]} = $line_AC;
			 }
			 else
			 {# this AC has more than one AC line, exists in multi-places
				 #print "This AC $v_AC[$i] has more than one AC lines, $AC2ACline{$v_AC[$i]} and  $line_AC  \n";
                 $AC2ACline{$v_AC[$i]}.="|".$line_AC;
			 }

		 }
		 if($v_gene eq "NULL")
		 {
			 #print "No gene for this protein : $v_ID , @v_AC\n";
		 }
		 else
		 {# we get geneID for this protein, now build relationship of uniprot ID and geneID
			 if(not exists $ID2gene{$v_ID})
			 {
				 $ID2gene{$v_ID}=$v_gene;
			 }
			 else
			 {
				 #print "this ID $v_ID has more than one geneIDs , one of them is $v_gene, another is $ID2gene{$v_ID}!\n";
			 }
			 if(not exists $gene2ID{$v_gene})
			 {
				 $gene2ID{$v_gene}=$v_ID;
			 }
			 else
			 {
				 #print "this gene ID $v_gene has more than one ID , one of them is $v_ID, another is $gene2ID{$v_gene}!\n";
			 }
			 if(not exists $gene2ACline{$v_gene})
			 {
				 $gene2ACline{$v_gene}=$line_AC;
			 }
			 else
			 {
				 #print "This gene ID $v_gene has more than one AC lines!\n";
                 $gene2ACline{$v_gene}.="|".$line_AC;
			 }
             
		 }

         $v_ID="NULL";
		 @v_AC=();
		 $index_AC=0;
		 $v_gene="NULL";
		 $line_AC="NULL";

	     next;
	 }
	 if($tem_split[0] eq "ID")
	 {
		 $v_ID=$tem_split[1];
	 }
	 elsif($tem_split[0] eq "AC")
	 {
		 for($i=1;$i<@tem_split;$i++)
		 {
			 @tem_2=split(/\;/,$tem_split[$i]);
			 $v_AC[$index_AC]=$tem_2[0];
			 $index_AC++;
		 }
		 $line_AC=$line;
	 }
	 elsif(($tem_split[0] eq "DR") && ($tem_split[1] eq "GeneID;"))
	 {# get the geneID
		 @tem_2=split(/\;/,$tem_split[2]);
		 $v_gene="GeneID:".$tem_2[0];
	 }

 }

 $IN->close();
 print "Done!\n";
 print "Loading protein lists ...";
 #############################################################################################
 #
 #                2. Load protein lists
 #
 #############################################################################################
 my(@protein_list)=();           # the list for all protein which needs to make prediction
 my(@protein_value)=();          # the evalue and target name for this protein
 my(@protein_type)=();           # the network type for this protein
 my(@protein_mapping)=();        # the mapped protein/gene for this protein 
 my($protein_index)=0;
 $IN=new FileHandle "$addr_list";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<1)
	 {
		 next;
	 }
	 $protein_list[$protein_index]=$tem_split[0];
	 $protein_value[$protein_index]=$tem_split[1];
	 $protein_type[$protein_index]="NULL";
	 $protein_mapping[$protein_index]="NULL";
	 $protein_index++;

 }
 $IN->close();
 print "Done, $protein_index proteins loaded!\n";
 print "Assigning network type for each protein ...";
 
 #############################################################################################
 #
 #                3. Load each network, decide the network type and mapped protein/geneID
 #
 #############################################################################################
 my(%hash_all_AClines)=();                # store all AC lines for this network, the key is the AC lines, and the value is the proteinID in the network, could be multiple, ID1|ID2 ...

 my(%unique_ID)=();                  #  store the unique ID for this network, in case we process the same protein/gene more than 1 times
 opendir(DIR,"$dir_networks"); 
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 $path_input=$dir_networks."/".$file;
	 @tem_split=split(/\./,$file);
	 if(@tem_split<2 || $tem_split[1] ne "network")
	 {
		 print "We get $file, but the name is not *.network, we skip this file!\n";
	 }
	 $network_type=$tem_split[0];

	 print "Processing $file ...";
	 %hash_all_AClines=(); 
     %unique_ID=();
	 

	 $IN=new FileHandle "$path_input";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split == 4)
		 {# this is protein protein interaction network, protein1 protein2 network_type network_type
 			 if(not exists $AC2ACline{$tem_split[0]})
			 {# this AC is a wrong AC ,could not found in the uniprot, check the reason
				 print "This AC $tem_split[0] is not a common AC number, not exists in the uniprot database, check this !\n";
				 #sleep(1);
				 next;
			 }
			 else
			 {
				 ###### check whether we process this before ########
				 if(exists $unique_ID{$tem_split[0]})
				 {
					 next;
				 }
				 else
				 {
					 $unique_ID{$tem_split[0]}=1;
				 }
				 $key=$AC2ACline{$tem_split[0]};
				 @tem_2=split(/\|/,$key);         # see whether this AC has multi AC lines
				 for($i=0;$i<@tem_2;$i++)
				 {
					 if(not exists $hash_all_AClines{$tem_2[$i]})
					 {
						 $hash_all_AClines{$tem_2[$i]}=$tem_split[0];  
					 }
					 else
					 {
						 $hash_all_AClines{$tem_2[$i]}.="|".$tem_split[0];
					 }
				 }
                 
			 }
 			 if(not exists $AC2ACline{$tem_split[1]})
			 {# this AC is a wrong AC ,could not found in the uniprot, check the reason
				 print "This AC $tem_split[1] is not a common AC number, not exists in the uniprot database, check this !\n";
				 #sleep(1);
				 next;
			 }
			 else
			 {
				 ###### check whether we process this before ########
				 if(exists $unique_ID{$tem_split[1]})
				 {
					 next;
				 }
				 else
				 {
					 $unique_ID{$tem_split[1]}=1;
				 }
				 $key=$AC2ACline{$tem_split[1]};
				 @tem_2=split(/\|/,$key);         # see whether this AC has multi AC lines
				 for($i=0;$i<@tem_2;$i++)
				 {
					 if(not exists $hash_all_AClines{$tem_2[$i]})
					 {
						 $hash_all_AClines{$tem_2[$i]}=$tem_split[1];  
					 }
					 else
					 {
						 $hash_all_AClines{$tem_2[$i]}.="|".$tem_split[1];
					 }
				 }
			 }
		 }
		 elsif(@tem_split==3)
		 {# this is a gene-gene interaction network. Gene1 number_of_interaction gene2
 			 if(not exists $gene2ACline{$tem_split[0]})
			 {# this AC is a wrong AC ,could not found in the uniprot, check the reason
				 print "This AC $tem_split[0] is not a common AC number, not exists in the uniprot database, check this !\n";
				 #sleep(1);
				 next;
			 }
			 else
			 {
				 ###### check whether we process this before ########
				 if(exists $unique_ID{$tem_split[0]})
				 {
					 next;
				 }
				 else
				 {
					 $unique_ID{$tem_split[0]}=1;
				 }
				 $key=$gene2ACline{$tem_split[0]};
				 @tem_2=split(/\|/,$key);         # see whether this AC has multi AC lines
				 for($i=0;$i<@tem_2;$i++)
				 {
					 if(not exists $hash_all_AClines{$tem_2[$i]})
					 {
						 $hash_all_AClines{$tem_2[$i]}=$tem_split[0];  
					 }
					 else
					 {
						 $hash_all_AClines{$tem_2[$i]}.="|".$tem_split[0];
					 }
				 }
                 
			 }
 			 if(not exists $gene2ACline{$tem_split[2]})
			 {# this AC is a wrong AC ,could not found in the uniprot, check the reason
				 print "This AC $tem_split[2] is not a common AC number, not exists in the uniprot database, check this !\n";
				 #sleep(1);
				 next;
			 }
			 else
			 {
				 ###### check whether we process this before ########
				 if(exists $unique_ID{$tem_split[2]})
				 {
					 next;
				 }
				 else
				 {
					 $unique_ID{$tem_split[2]}=1;
				 }
				 $key=$AC2ACline{$tem_split[2]};
				 @tem_2=split(/\|/,$key);         # see whether this AC has multi AC lines
				 for($i=0;$i<@tem_2;$i++)
				 {
					 if(not exists $hash_all_AClines{$tem_2[$i]})
					 {
						 $hash_all_AClines{$tem_2[$i]}=$tem_split[2];  
					 }
					 else
					 {
						 $hash_all_AClines{$tem_2[$i]}.="|".$tem_split[2];
					 }
				 }
			 }		 
		 }
		 else
		 {
			 print "Warning, we skip $line! We don't know the interaction type for this network!\n";
			 next;
		 }
	 }
	 $IN->close();
     ###### now we already get all AC lines for this network, check the protein list #######
     for($i=0;$i<$protein_index;$i++)
	 {# fill in the protein list, type , mapping
		 if (not exists $AC2ACline{$protein_list[$i]})
		 {# this protein is not in the uniprot ???
			 print "This protein ID $protein_list[$i] is not in the uniprot , check the reason!\n";
			 next;
		 }
		 $key=$AC2ACline{$protein_list[$i]};          # get the AC line
		 if(exists $hash_all_AClines{$key})
		 {# this protein is mapped to the processing network, good
			 if($protein_type[$i] eq "NULL")
			 {# first time assign the network type 
				 $protein_type[$i] = $network_type;           # file name is the network type
			 }
			 else
			 {# this protein has mapping to multi networks
				 $protein_type[$i].=";".$network_type;
			 }
             #@tem_split=split(/\|/,$hash_all_AClines{$key});  # this is the mapped protein ID in the network
			 if($protein_mapping[$i] eq "NULL")
			 {
				 $protein_mapping[$i]=$hash_all_AClines{$key};
			 }
			 else
			 {# this protein map to multi network
				 $protein_mapping[$i].=";".$hash_all_AClines{$key};
			 }
		 }
	 }
 }
 $key=0;
 for($i=0;$i<$protein_index;$i++)
 {
	 if($protein_type[$i] eq "NULL")
	 {
		 $key++;
	 }
 }
 print "Done! $key proteins out of total $protein_index are not mapped to any network!\n";

 #############################################################################################
 #
 #                4. Output the result
 #    Format is : 
 #    Column 1: proteinID
 #    Column 2: evalue and target name , such as 1e-10|T01234
 #    Column 3: network type mapped, could map to multi networks, such as 562;3_contact_NORMAL_B_gene;4932
 #    Column 4: mapped protein/geneID in the network, could have multiID in multi-networks. Such as P0CY12|P0CY10;Q9VR90|Q9VHG6|Q24297;C6TP87
 #
 #############################################################################################  
 $OUT = new FileHandle ">$addr_output";
 for($i=0;$i<$protein_index;$i++)
 {
	 print $OUT $protein_list[$i]."\t".$protein_value[$i]."\t".$protein_type[$i]."\t".$protein_mapping[$i]."\n";

 }
 
 $OUT->close();

  my $endtime = localtime();
  print "\nThe time ended at : $endtime.\n";

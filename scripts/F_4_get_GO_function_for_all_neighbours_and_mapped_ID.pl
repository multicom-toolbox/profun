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



  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
      print("This program load the protein lists with network type and neighbors, search uniprot database for the GO function of all proteins/genes ...\n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/7/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_uniprot_dat protein_list(has neighbor information) addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../result/F_3_protein_lists_with_network_and_one_step_neighbour ../result/F_4_GO_function_for_all_protein_gene\n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";
 
 my($addr_uniprot)=$ARGV[0];
 my($addr_list)=$ARGV[1];      # protein list
 my($output)=$ARGV[2];         # get the function for all protein or gene used
 $| = 1;
 print "Loading protein lists ...";
##############################################################################
#
#
#               1. Load all protein ID or gene ID
#
#
##############################################################################
  my(%hash_all_IDs)=();
  my($IN,$OUT,$line,$i,$j,$key);
  my(@tem_split,@tem_2,@tem_3);
  $IN=new FileHandle "$addr_list";
  while(defined($line=<$IN>))
  {
     chomp($line);
     @tem_split=split(/\s+/,$line);
     if(@tem_split<5)
     {
        print "Warning! The input file should have 5 columns, proteinID evalue_target_name species mapped_protein/gene protein/gene_neighbours\n";
        next;
     }
     if(not exists $hash_all_IDs{$tem_split[0]})     # get the query proteinID
     {
        $hash_all_IDs{$tem_split[0]}="NULL";
     }
     @tem_2=split(/\;/,$tem_split[3]);       # get the mapped proteinID/geneID
     for($i=0;$i<@tem_2;$i++)
     {
         if(not exists $hash_all_IDs{$tem_split[0]})
         {
             $hash_all_IDs{$tem_split[0]}="NULL";
         }
     }
     @tem_2=split(/\|/,$tem_split[4]);      # get the neighbors of mapped protein/gene
     for($i=0;$i<@tem_2;$i++)
     {# get the neighbors for mapped protein i
         @tem_3=split(/\_/,$tem_2[$i]);
         for($j=0;$j<@tem_3;$j++)
         {
             if(not exists $hash_all_IDs{$tem_3[$j]})
             {
                 $hash_all_IDs{$tem_3[$j]}="NULL";
             }
         }
     }
     
  }
  $IN->close(); 
  print "We get all protein/gene IDs to search against uniprot database ...\n";

##############################################################################
#
#
#               2. Search uniprot to get the GO function for each protein/ gene
#
#
##############################################################################
  my($geneID,$GO_vector,$AC_line);
  
  
  $geneID="NULL";           # store the geneID for this record
  $AC_line="NULL";          # store all the AC numbers
  $GO_vector="NULL";        # store the GO functions
  $IN=new FileHandle "$addr_uniprot";
  while(defined($line=<$IN>))
  {
	  chomp($line);
	  $line=~s/\s+$//;
	  @tem_split=split(/\s+/,$line);
	  if($line eq "//")
	  {# we come to the end of this record
		  if($GO_vector ne "NULL")
		  {
			  if(exists $hash_all_IDs{$geneID})
			  {# add the geneID 's function
				  $hash_all_IDs{$geneID}=$GO_vector;
			  }
			  @tem_2=split(/\s+/,$AC_line);
              for($i=1;$i<@tem_2;$i++)
			  {
				  @tem_3=split(/\;/,$tem_2[$i]);
				  if(exists $hash_all_IDs{$tem_3[0]})
				  {# add each AC protein
					  $hash_all_IDs{$tem_3[0]}=$GO_vector;
				  }
			  }
		  }
          $GO_vector="NULL";
		  $AC_line="NULL";
		  $geneID="NULL";
	  }
	  if($tem_split[0] eq "AC")
	  {# this is AC line
		  $AC_line=$line;
	  }
	  elsif($tem_split[0] eq "DR")
	  {
		  if($tem_split[1] eq "GeneID;")
		  {# this is a geneID
			  @tem_2=split(/\;/,$tem_split[2]);
			  $geneID="GeneID:".$tem_2[0];
		  }
		  elsif($tem_split[1] eq "GO;")
		  {
			  @tem_2=split(/\;/,$tem_split[2]);

			  if($GO_vector eq "NULL")
			  {
				  $GO_vector=$tem_2[0];
			  }
			  else
			  {
				  $GO_vector.="_".$tem_2[0];
			  }
		  }
	  }
  }
  $IN->close();
##############################################################################
#
#
#               3. Output the GO function for protein/gene
#
#
##############################################################################
  $OUT=new FileHandle ">$output";
  foreach $key (keys %hash_all_IDs)
  {
	  if($hash_all_IDs{$key} ne "NULL")
	  {# this protein/gene has function
		  print $OUT $key."\t".$hash_all_IDs{$key}."\n";
	  }
  }
  $OUT->close();
#We search uniprot database to get functions for all one step neighbours( convert geneID or protein AC to AC lines, and search GO for that AC lines). And also search functions for the mapped protein/gene, since we have another method which uses the multiple hits statistis for making function prediction.


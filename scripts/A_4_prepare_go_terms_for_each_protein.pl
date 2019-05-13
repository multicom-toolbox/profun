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

      print "This program use the ppi network as input, and extact the go term from the uniprot database to generates the go terms for each protein!\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 6/27/2013\n\n");
	  print("***********************************************************\n");
	  print "perl $0 gene_network uniprot_database output\n";
	  print "For example:\n";
	  print "perl $0 ../data/converted_network_with_taxid/other.network /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../result/2_real_go_function_for_each_protein/other.mitab.go_functions\n";


	  exit(0);
 }
 my($input_network)=$ARGV[0];
 my($input_uniprot)=$ARGV[1];
 my($output)=$ARGV[2];
 -s $input_network || die "cannot open input gene gene network $input_network\n";
 -s $input_uniprot || die "cannot open input uniprot database $input_uniprot\n";
 
 my($key,$value,$IN,$OUT,$line,$i,$j);
 my(@tem_split,@tem_2);
 my(%hash_all_geneID)=();                   # store all geneID used in the highly contact gene pairs

 #################################################################
 #
 #           1. Get all protein ID
 #
 #################################################################
 $IN=new FileHandle "$input_network";
 defined($IN) || die "Cannot open input $input_network!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<4)
	 {# this is a ppi network file, so has protein1 protein2 taxid1 taxid2
		 next;
	 }
	 if(!exists $hash_all_geneID{$tem_split[0]})
	 {
		 $hash_all_geneID{$tem_split[0]}=1;
	 }
     if(!exists $hash_all_geneID{$tem_split[1]})
	 {
		 $hash_all_geneID{$tem_split[1]}=1;
	 }
 }
 $IN->close();
 print "All protein ID has been loaded!\n";

 #################################################################
 #
 #           2. Search swiss prot to get all GO functions for all genes
 #
 ################################################################# 
    my(%gene_GO)=();
    
	my ($to_search)=1; # 1 means it need to search
	my ($hit_gene) = 0; # to see whether this protein is what we want
	my($hit_go)=0;
    
	my($hit_seq)=0;
	my($seq_vector)="NULL";

	my($GO_vector)=":";
	my($saved_gene);
    $IN = new FileHandle "$input_uniprot";
    if (! defined($IN)) 
    {
        print("Can't open spec file $input_uniprot: $!\n");
       return;
    }
	while ( defined($line = <$IN>))
	{
	   $line =~ s/\s+$//;
       chomp($line); # strip off trailing '\n' (newline)
       if($line eq "//")
	   {
############### save the GO function for each AC#####################
           if($hit_gene==1 && $hit_go==1 )
		   {
			   if(! exists $gene_GO{$saved_gene})
			   {
				   $gene_GO{$saved_gene}=$GO_vector;
			   }
			   ########## save sequence for each AC #########
#			   if(! exists $hash_sequence_of_geneID{$saved_gene})
#			   {
#				   print "The gene is not in the hash list???? Error!\n";
#			   }
#			   $hash_sequence_of_geneID{$saved_gene}=$seq_vector;
		   }
		   



############### save the GO function for each AC#####################
           $hit_seq=0;
		   $seq_vector="NULL"; 

		   $to_search = 1;
		   $hit_gene=0;
		   $hit_go=0;
		   $GO_vector=":";
		   next;
	   }
	   if($to_search == 1)
	   {#we will search
		   @tem_split = split(/\s+/,$line);
		   if($tem_split[1] eq "GO;")
		   {
			   @tem_split=split(/\;/,$line);
			   if($GO_vector eq ":")
			   {
				   $i=1;
				   @tem_2=split(/\s+/,$tem_split[$i]);
				   $GO_vector=$tem_2[1];
				   for($i=4;$i<@tem_split;$i+=3)
				   {
					   @tem_2=split(/\s+/,$tem_split[$i]);
                       $GO_vector=$GO_vector."|".$tem_2[1];
				   }
			   }
			   else
			   {
				   for($i=1;$i<@tem_split;$i+=3)
				   {
					   @tem_2=split(/\s+/,$tem_split[$i]);
                       $GO_vector=$GO_vector."|".$tem_2[1];
				   }
			   }
               $hit_go=1;
		   }
		   elsif($tem_split[0] eq "AC")
		   {  # we get the protein AC
                        for($i=1;$i<@tem_split;$i++)
                        {
			   @tem_2=split(/\;/,$tem_split[$i]);
			   $saved_gene=$tem_2[0];
			   if(exists $hash_all_geneID{$saved_gene})
			   {
                                $hit_gene=1;
                                last;
			   }
                        }
			   if($i==@tem_split)
			   {# we do not hit a gene
				   $to_search=0;
			   }
		   }
		   elsif($tem_split[0] eq "SQ")
		   {# this is for getting the sequence
			   $hit_seq=1;
		   }
		   else
		   {
			   if($hit_seq == 1)
			   {# this is for sequence 
				   if($seq_vector eq "NULL")
				   {
#					   $seq_vector=add_sequence($line);
				   }
				   else
				   {
#					   $seq_vector.=add_sequence($line);      # add the sequence to the sequence vector, if seq_vector is "NULL", then initialize it, otherwise, add the new line.

				   }
				   
			   }
		   }
	   }
	}#end while
    $IN->close();
print "Swiss prot database checked!\n";
  ########## output the gene function ###############
  ########## we only keep the GO term GO:000***, since the fastsemsim only accept this ########
 
  $OUT = new FileHandle ">$output";
  my($has_go);

  foreach $key (keys %gene_GO) 
  {
	  @tem_split=split(/\|/,$gene_GO{$key});
	  $has_go=0;              # tag for whether the gene has go terms
	  for($i=0;$i<@tem_split;$i++)
	  {
		  if($has_go == 0)
		  {
			  print $OUT $key."|".$tem_split[$i];
			  $has_go=1;
		  }
		  else
		  {
			  print $OUT "\t".$tem_split[$i];
		  }
	  }
	  if($has_go == 1)
	  {
		  print $OUT "\n";
	  }
  }
  $OUT->close();

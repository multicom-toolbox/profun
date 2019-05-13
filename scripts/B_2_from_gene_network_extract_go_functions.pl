##########################################################################################################
#    Function about getting all GO terms of contact gene pairs,                     check GO in clster.  #
#                                  output: addr_output                                                   #
#										Renzhi Cao  													 #            
#																										 #            
#									    11/25/2012														 #            
#																										 #            
#																										 #            
#									Revised at 4/12/2013												 #            
#																										 #            
##########################################################################################################
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
 sub add_sequence($);             # add a sequence to a sequence vector.

  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args
          print "Revised at 10/3/2013, we add the number of interactions of that two genes to each GO pairs!\n";
	  print "Revised at 4/12/2013, we add the geneID and gene sequence to the sentence #ENDOFGENE\t#ENDOFGENE here \n";
	 print "revised at 7/30/2013, add a function, if number of selected gene pairs is less than 0, than we select all gene pairs!\n";
           print("This program will select gene pairs with a threshold contact, also check uniprot to make sure all gene has GO functions. output the GO functions of these gene pairs. Output format is GO:***** GO:****** for each line, and #ENDOFGENE,#ENDOFCLUSTER.\n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 11/25/2012\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_swiss_prot_dat addr_input_gene-gene-contact-network(The file has format like:geneID1 5 geneID2 for each row, the contact between each gene pair is more than a threshold) number_of_gene_pairs_to_select(negative value to select all gene pairs) addr_output!\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../data/gene_interaction_network/3_contact_gene_pairs_NORMAL_B -1 ../result/B_2_gene_network_go_function_pairs/all_go_pairs_for_3_NORMAL\n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";
 my($addr_swiss)=$ARGV[0];
 my($input)=$ARGV[1];
 my($total_num)=$ARGV[2];
 my($addr_output)=$ARGV[3];

 -s $input || die "No input file $input!\n";
 my($key,$value,$IN,$OUT,$line,$i,$j);
 my(@tem_split,@tem_2);
 my(%hash_all_geneID)=();                   # store all geneID used in the highly contact gene pairs

 my(%hash_sequence_of_geneID)=();           # store all sequence of the geneID

 my(%hash_interaction)=();                # store the number of interactions between two genes
 


 #################################################################
 #
 #           1. Get all geneID
 #
 #################################################################
 $IN=new FileHandle "$input";
 defined($IN) || die "Cannot open input $input!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {# this is a gene network file, so has gene1 number gene2
		 next;
	 }
	 if(!exists $hash_all_geneID{$tem_split[0]})
	 {
		 $hash_all_geneID{$tem_split[0]}=1;
	 }
     if(!exists $hash_all_geneID{$tem_split[2]})
	 {
		 $hash_all_geneID{$tem_split[2]}=1;
	 }

	 if(!exists $hash_sequence_of_geneID{$tem_split[0]})
	 {
		 $hash_sequence_of_geneID{$tem_split[0]}="NULL";
	 }
     if(!exists $hash_sequence_of_geneID{$tem_split[2]})
	 {
		 $hash_sequence_of_geneID{$tem_split[2]}="NULL";
	 }


        $key=$tem_split[0]."_".$tem_split[2];
        if(not exists $hash_interaction{$key})
        {
            $hash_interaction{$key}=$tem_split[1];
        }
        else 
        {next;}
        $key=$tem_split[2]."_".$tem_split[0];
        if(not exists $hash_interaction{$key})
        {
            $hash_interaction{$key}=$tem_split[1];
        }
 }
 $IN->close();
 print "All gene ID has been loaded!\n";
 #################################################################
 #
 #           2. Search swiss prot to get all GO functions for all genes
 #
 ################################################################# 
    my(%gene_GO)=();
    my(@tem_2);
	my ($to_search)=1; # 1 means it need to search
	my ($hit_gene) = 0; # to see whether this protein has the gene we want
	my($hit_go)=0;
	my($hit_OS)=0;
    
	my($hit_seq)=0;
	my($seq_vector)="NULL";

	my($GO_vector)=":";
	my($saved_gene);
    $IN = new FileHandle "$addr_swiss";
    if (! defined($IN)) 
    {
        print("Can't open spec file $addr_swiss: $!\n");
       return;
    }
	while ( defined($line = <$IN>))
	{
	   $line =~ s/\s+$//;
       chomp($line); # strip off trailing '\n' (newline)
       if($line eq "//")
	   {
############### save the GO function for each AC#####################
           if($hit_gene==1 && $hit_go==1 && $hit_OS ==1)
		   {
			   if(! exists $gene_GO{$saved_gene})
			   {
				   $gene_GO{$saved_gene}=$GO_vector;
			   }
			   ########## save sequence for each AC #########
			   if(! exists $hash_sequence_of_geneID{$saved_gene})
			   {
				   print "The gene is not in the hash list???? Error!\n";
			   }
			   $hash_sequence_of_geneID{$saved_gene}=$seq_vector;
		   }
		   



############### save the GO function for each AC#####################
           $hit_seq=0;
		   $seq_vector="NULL"; 

		   $to_search = 1;
		   $hit_gene=0;
		   $hit_go=0;
		   $hit_OS=0;
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
		   elsif($tem_split[0] eq "OS")
		   {
               if($tem_split[2] eq "sapiens")
			   {
				   $hit_OS=1;
			   }
		       if($hit_OS==0)
			   {# we do not hit a OS
				   $to_search=0;
			   }
		   }
		   elsif($tem_split[1] eq "GeneID;")
		   {  
			   @tem_2=split(/\;/,$tem_split[2]);
			   $saved_gene="GeneID:".$tem_2[0];
			   if(exists $hash_all_geneID{$saved_gene})
			   {
                   $hit_gene=1;
			   }
			   else
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
					   $seq_vector=add_sequence($line);
				   }
				   else
				   {
					   $seq_vector.=add_sequence($line);      # add the sequence to the sequence vector, if seq_vector is "NULL", then initialize it, otherwise, add the new line.

				   }
				   
			   }
		   }
	   }
	}#end while
    $IN->close();
print "Swiss prot database checked!\n";
 #################################################################
 #
 #           3. Read gene pairs file again, randomly select $total_num of gene pairs which have function,output the GO term
 #
 #################################################################
  if(-e $addr_output)
  {
    print "the result file | : $addr_output  ...Exists!\n"; 
  }
  else 
  { 
    open (File, "&gt;$addr_output");
    chmod (0777, $addr_output); 
    close (File);
  }      
  $OUT = new FileHandle "> $addr_output";
  defined($OUT) || die "Cannot open output file $addr_output!\n";

 my(@array_contact_pairs)=();            # store gene pairs
 my(@go_1,@go_2,@filtered_go_1,@filtered_go_2);
 my($index_array)=0;
 $IN=new FileHandle "$input";
 defined($IN) || die "Cannot open input $input!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {# this is a gene network file, so has gene1 number gene2
		 next;
	 }
     $array_contact_pairs[$index_array]=$tem_split[0]."|".$tem_split[2];
	 $index_array++;
 }
 $IN->close();
 if($total_num>=$index_array)
 {# select total_num, larger equal to all gene pairs, so directly select all of them
	 print "We only have $index_array gene pairs, so set $total_num to $index_array!\n";
	 $total_num=$index_array;
 }
 if($total_num<0)
 { # add this, so we will select all gene pairs.
      $total_num=$index_array; 
 }
 my($key_num_interaction);
 my($random_range)=$index_array;
 my($random_total)=0;              # count the random number of gene pairs we select
 my($random_selected,$random_selected_pair,$genei,$genej,$k1,$k2);
 while($random_total<$total_num)
 {#
	 $random_selected=int(rand($random_range));
     $random_selected_pair=$array_contact_pairs[$random_selected];    # the gene pair selected
     for($i=$random_selected;$i<$random_range-1;$i++)
	 {
        $array_contact_pairs[$i]=$array_contact_pairs[$i+1];
	 }
     $random_range--;
     if($random_range<0)
	 {
		print "Not enough pair of genes have GO terms for you. Now get $random_total!\n";
		last;
	 }
	 @tem_split=split(/\|/,$random_selected_pair);
	 $genei=$tem_split[0];
	 $genej=$tem_split[1];
         $key=$genei."_".$genej;
         if(not exists $hash_interaction{$key})
         {
             print "This should not happen, we get two genes : $key, and they are actually not in contact, why???\n";
             exit(0);
         }
         $key_num_interaction=$hash_interaction{$key};
     if((exists $gene_GO{$genei}) && (exists $gene_GO{$genej}))
	 {# this gene pair has GO function
		 # output the GO function for this two genes
         @go_1=split(/\|/,$gene_GO{$genei});
         @go_2=split(/\|/,$gene_GO{$genej});
		 ########### only keep the GO:000***, because the tool FastSemsim only accept this #########
         @filtered_go_1=();
		 @filtered_go_2=();
		 $k2=0;
		 for($k1=0;$k1<@go_1;$k1++)
		 {
			 if(substr($go_1[$k1],0,6) eq "GO:000")
			 {
				 $filtered_go_1[$k2]=$go_1[$k1];
				 $k2++;
			 }
		 }
		 if($k2==0)
		 {# no function we need
			 next;
		 }
		 $k2=0;
		 for($k1=0;$k1<@go_2;$k1++)
		 {
			 if(substr($go_2[$k1],0,6) eq "GO:000")
			 {
				 $filtered_go_2[$k2]=$go_2[$k1];
				 $k2++;
			 }
		 }
		 if($k2==0)
		 {# no function we need
			 next;
		 }
		 ###########################################################################################
         ##### output the pairwise GO function ######
		 for($k1=0;$k1<@filtered_go_1;$k1++)
		 {
			 for($k2=0;$k2<@filtered_go_2;$k2++)
			 {
				 print $OUT $filtered_go_1[$k1]."\t".$filtered_go_2[$k2]."\t".$key_num_interaction."\n";
			 }
		 }
         
         print $OUT "#ENDOFGENE\t#ENDOFGENE $genei|$hash_sequence_of_geneID{$genei}|$genej|$hash_sequence_of_geneID{$genej}\n";
		 $random_total++;
	 }
	 else
	 {# at least one of them has no GO function, skip
		 next;
	 }
 }
 $OUT->close();



 sub add_sequence($)
 {# extract the sequence from a line. 
	 my($line)=@_;
	 my(@tem_split)=split(/\s+/,$line);
	 if($tem_split[0] ne "")
	 {
		 print "Check this line : $line\n";
		 print "Maybe the sequence is wrong??? \n";
		 exit(0);
	 }
	 my($seq);
	 $seq=$tem_split[1];
	 my($i);
	 for($i=2;$i<@tem_split;$i++)
	 {
		 $seq.=$tem_split[$i];
	 }
	 $seq=uc($seq);                # convert to uppercase
	 return $seq;
 }

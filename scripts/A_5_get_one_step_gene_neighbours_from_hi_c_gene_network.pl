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
 if (@ARGV != 2)
 { # @ARGV used in scalar context = number of args

      print "This program use the gene network as input, get one step gene neighbours for each gene, we don't consider the gene contact with itself, should we? !\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 6/27/2013\n\n");
	  print "Revised by Renzhi at 7/1/2013, the gene neighbours are ranked by the number of contact, so we can use this information to make prediction\n";
	  print("***********************************************************\n");
	  print "perl $0 gene_network output_gene_neighbours\n";
	  print "For example:\n";
          print "perl $0 ../data/gene_interaction_network/3_contact_gene_pairs_NORMAL_B ../result/3_one_step_protein_neighbours/3_contact_NORMAL_B_gene.one_step_gene_neighbours\n";
          exit(0);
 }
 my($input_network)=$ARGV[0];
 my($output)=$ARGV[1];

 my($IN,$line,$OUT,$key,$neighbour,$i,$j,$for_rank,$c1,$c2);
 my(@tem_split,@neighbours);
 my(@contacts);


 my(%hash_contact)=();        # store the geneID1|geneID2 and the number of contact between them

 ######### get all gene neighbours for each gene ######
 my(%hash)=();

 $IN=new FileHandle "$input_network";
 defined($IN) || die "Cannot open input $input_network!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {# this is a gene network file, so has gene1 number gene2
		 next;
	 }
	 if($tem_split[0] eq $tem_split[2])
	 {# we don't use the gene itself as the neighbour
		 next;
	 }
	 $key=$tem_split[0]."|".$tem_split[2];
	 if(not exists $hash_contact{$key})
	 {# check whether we need to add this gene pairs and the contact
         $hash_contact{$key}=$tem_split[1];
		 $key=$tem_split[2]."|".$tem_split[0];
         $hash_contact{$key}=$tem_split[1]; 
	 }
     
	 if(not exists $hash{$tem_split[0]})
	 {# we get a neighbour for $tem_split[0]
		 if($tem_split[0] ne $tem_split[2])
		 {# skip the gene itself
			 $hash{$tem_split[0]}=$tem_split[2];
		 }
	 }
	 else
	 {# add the new neighbour
		 ###### insert the new neighbour, and keep the highest contact neighbour rank at first #######
         ### first check whether it already has this neighbour #####
		 @neighbours = split(/\s+/,$hash{$tem_split[0]});
         for($i=0;$i<@neighbours;$i++)
		 {
			 if($neighbours[$i] eq $tem_split[2])
			 {# we already add this gene as the neighbour
				 last;
			 }
		 }
		 $neighbours[$i]=$tem_split[2];     # add this new neighbour as the neighbour of tem_split[0]
		 ###### now use insert ranking to rank the array #######
		 $for_rank=$neighbours[@neighbours-1];      # get the last element, insert this again to rank it
		 $j=@neighbours-1;                            # get the last element index
		 for(;$j>0;$j--)
		 {
			 $key=$tem_split[0]."|".$neighbours[$j];
			 if(not exists $hash_contact{$key})
			 {
				 die "Check the $key is not exists for the number of contact in the gene network, why 1???\n";
			 }
			 $c1=$hash_contact{$key};
			 $key=$tem_split[0]."|".$neighbours[$j-1];
			 if(not exists $hash_contact{$key})
			 {
				 die "Check the $key is not exists for the number of contact in the gene network, why 11???\n";
			 }
			 $c2=$hash_contact{$key};
			 if($c2>=$c1)
			 {
				 last;
			 }
			 else
			 {
				 $neighbours[$j]=$neighbours[$j-1];
			 }
		 }
		 $neighbours[$j]=$for_rank;
         
		 $hash{$tem_split[0]}="@neighbours";
	 }


	 if(not exists $hash{$tem_split[2]})
	 {# we get a neighbour for $tem_split[2]
		 if($tem_split[2] ne $tem_split[0])
		 {# skip the gene itself
			 $hash{$tem_split[2]}=$tem_split[0];
		 }
	 }
	 else
	 {# add the new neighbour
		 ###### insert the new neighbour, and keep the highest contact neighbour rank at first #######
         ### first check whether it already has this neighbour #####
		 @neighbours = split(/\s+/,$hash{$tem_split[2]});
         for($i=0;$i<@neighbours;$i++)
		 {
			 if($neighbours[$i] eq $tem_split[0])
			 {# we already add this gene as the neighbour
				 last;
			 }
		 }
		 $neighbours[$i]=$tem_split[0];     # add this new neighbour as the neighbour of tem_split[0]
#print "$tem_split[2]'s neighbours: @neighbours\t";
		 ###### now use insert ranking to rank the array #######
		 $for_rank=$neighbours[@neighbours-1];      # get the last element, insert this again to rank it
		 $j=@neighbours-1;                            # get the last element index
		 for(;$j>0;$j--)
		 {
			 $key=$tem_split[2]."|".$neighbours[$j];
			 if(not exists $hash_contact{$key})
			 {
				 die "Check the $key is not exists for the number of contact in the gene network, why 2 ???\n";
			 }
			 $c1=$hash_contact{$key};
			 $key=$tem_split[2]."|".$neighbours[$j-1];
			 if(not exists $hash_contact{$key})
			 {
				 die "Check the $key is not exists for the number of contact in the gene network, why  22???\n";
			 }
			 $c2=$hash_contact{$key};
			 if($c2>=$c1)
			 {
				 last;
			 }
			 else
			 {
				 $neighbours[$j]=$neighbours[$j-1];
			 }
		 }
		 $neighbours[$j]=$for_rank;
#print "@neighbours\n";
		 $hash{$tem_split[2]}="@neighbours";
	 }

 }
 $IN->close();
 #### output the gene neighbours ####
 $OUT = new FileHandle ">$output";
 foreach $key (keys %hash) 
 {
	 print $OUT $key."|".$hash{$key}."\n";
 }
 $OUT->close();


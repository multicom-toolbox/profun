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
 use List::Util qw/shuffle/;


 if (@ARGV != 5)
 { # @ARGV used in scalar context = number of args
      print "New, efficient!\n";
      print "This program makes protein function prediction from all protein neighbours which has go function, and use the go function of the neighbour to predict the go term of this protein. \n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 7/1/2013\n\n");
	  print("***********************************************************\n");
	  print "perl $0 go_pairs_statistic_frequency(from this get which go pairs usually exists together) gene_go_term(get the real go terms for each protein) one_step_protein_neighbours(get the neighbours for each protein) protein_list_for_making_prediction output_protein_function_prediction\n";
	  print "For example:\n";
	  print "perl $0 ../result/1_go_frequency/other.go_frequency ../result/2_real_go_function_for_each_protein/other.mitab.go_functions ../result/3_one_step_protein_neighbours/other.one_step_protein_neighbours ../result/2_real_go_function_for_each_protein/other.mitab.go_functions ../test/4_protein_function_predictions_from_one_step_neighbours/other.function_prediction\n";
          exit(0);
 }
 my($input_go_fre)=$ARGV[0];
 my($input_gene_go)=$ARGV[1];
 my($input_gene_neighbour)=$ARGV[2];
 my($input_gene_list)=$ARGV[3];          # gene list for making prediction
 my($output)=$ARGV[4];



 my($NUM_prediction)=10;             # the total number of predictions we need.

 -s $input_go_fre || die "cannot open $input_go_fre\n";
 -s $input_gene_go || die "cannot open $input_gene_go\n";
 -s $input_gene_neighbour || die "cannot open $input_gene_neighbour\n";
 -s $input_gene_list || die "cannot open $input_gene_list\n";

 my($IN,$OUT,$line);
 my(@tem_split);
 ############ 1. Load the go term statistical information ###############
 my(%go_fre)=();            # key is go term, value is go_term1_frequency1|go_term2_frequency2 ..., already ranked, we can see which go term exists more frequently
 $IN=new FileHandle "$input_go_fre";
 defined($IN) || die "Cannot open input $input_go_fre!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
     if(@tem_split<2)
	 {
		 print "skip $line\n";
		 next;
	 }
	 if(exists $go_fre{$tem_split[0]})
	 {
		 print "The go term $tem_split[0] is duplicated! ???\n";
	     exit(0);
	 }
	 $go_fre{$tem_split[0]}=$tem_split[1];
 }
 $IN->close();
print "Go statistics loaded!\n";
 ##########2. load the real go function for each protein ##############
 my(%gene_go)=();            # key is proteinID, value is go_term1 go_term2  ..., all the go terms for that proteinID
 $IN=new FileHandle "$input_gene_go";
 defined($IN) || die "Cannot open input $input_gene_go!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\|/,$line);
     if(@tem_split<2)
	 {
		 print "skip $line\n";
		 next;
	 }
	 if(exists $gene_go{$tem_split[0]})
	 {
		 print "The proteinID $tem_split[0] is duplicated! ???\n";
	     exit(0);
	 }
	 $gene_go{$tem_split[0]}=$tem_split[1];
 }
 $IN->close();
 print "Real go functions loaded!\n";
 ####### 3. Load the protein neighbours for each protein #################
 my(%gene_neighbour)=();            # key is proteinID, value is proteinID1 proteinID2 ..., all the neighbours of that proteinID
 $IN=new FileHandle "$input_gene_neighbour";
 defined($IN) || die "Cannot open input $input_gene_neighbour!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\|/,$line);
     if(@tem_split<2)
	 {
		 print "skip $line\n";
		 next;
	 }
	 if(exists $gene_neighbour{$tem_split[0]})
	 {
		 print "The proteinID $tem_split[0] is duplicated! ???\n";
	     exit(0);
	 }
	 $gene_neighbour{$tem_split[0]}=$tem_split[1];
 }
 $IN->close();
 print "Protein neighbours loaded!\n";
 ####### 4. Load the protein lists to make prediction #################
 my(@gene_lists)=();
 my($index)=0;
 $IN=new FileHandle "$input_gene_list";
 defined($IN) || die "Cannot open input $input_gene_list!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\|/,$line);
     if(@tem_split<1)
	 {
		 next;
	 }
	 $gene_lists[$index]=$tem_split[0];
	 $index++;
 }
 $IN->close();
 print "Protein lists loaded!\n";
 
 ##### for testing, we only make prediction for 1000 proteins ######
 @tem_split=split(/\//,$input_gene_list);
 if($tem_split[@tem_split-2] eq "2_real_go_function_for_each_protein")
 {
     print "Warning !!! We are predicting the protein function for testing, so only make prediction for 200 proteins! \n";
     if($index>200)
     {
         $index=200;
     }
 }
 
####### now make prediction for the input protein lists and output the result ########
 my($i,$j,$k,$l,$i2,$key);
 my(@all_neighbour_go)=();             # neighbours go term
 my(@all_neighbour_go_fre)=();         # frequency of neighbour go term

 my(@tem_neighbour);
 my($i_neighbour);
 my(%hash_GO_from_neighbour)=();
 my($index_all_neighbour)=0;
 $OUT = new FileHandle ">$output";

 my(@go_split1,@go_split2);
 for($i=0;$i<$index;$i++)
 {# make prediction for gene $gene_lists[$i] from it's neighbours go term
	 ### 1. Get one gene, with it's neighbour which has go term #######
	 if(not exists $gene_neighbour{$gene_lists[$i]})
	 {
		 print "Warning, this gene $gene_lists[$i] doesn't have a gene neighbour!!!\n";
		 next;
	 }
print "Predicting for $gene_lists[$i] ...\n";
	 ### 2. select all neighbours which has go term
	 @tem_neighbour=split(/\s+/,$gene_neighbour{$gene_lists[$i]});         # this contains all the gene neighbours
     %hash_GO_from_neighbour=();                                           # the GO from the neighbours

     for($i_neighbour=0;$i_neighbour<@tem_neighbour;$i_neighbour++)
	 {
	     ### 3. Already get one neighbour, use this neighbour's go to predict it's go term ####
		 if(not exists $gene_go{$tem_neighbour[$i_neighbour]})
		 {# this neighbour doesn't have GO function to use
			 next;
		 }
         @tem_split=split(/\s+/,$gene_go{$tem_neighbour[$i_neighbour]});      # this is all go terms of the neighbour
         for($j=0;$j<@tem_split;$j++)
		 {
			 ######## we get one GO of the neighbour, now find all partner of that GO and the frequency #########
			 if(exists $go_fre{$tem_split[$j]})
			 {# find the go term 
		         @go_split1=split(/\|/,$go_fre{$tem_split[$j]});
			     for($l=0;$l<@go_split1;$l++)
			     {
					 @go_split2=split(/\_/,$go_split1[$l]);    # get the go term and the frequency
                     if(not exists $hash_GO_from_neighbour{$go_split2[0]})
					 {
						 $hash_GO_from_neighbour{$go_split2[0]}=$go_split2[1];
					 }
					 else
					 {
						 $hash_GO_from_neighbour{$go_split2[0]}+=$go_split2[1];
					 }
				 }
			 }

		 }
	 }

    
     ########## rank and get the final output ###########
	 @go_split1=();
	 @go_split2=();
     $i_neighbour=0;
	 foreach $key (keys %hash_GO_from_neighbour)
	 {
		 for($i2=$i_neighbour;$i2>0;$i2--)
		 {
			 if($hash_GO_from_neighbour{$key}>$go_split2[$i2-1])
			 {
				 $go_split1[$i2]=$go_split1[$i2-1];
				 $go_split2[$i2]=$go_split2[$i2-1];
			 }
			 else
			 {
				 last;
			 }
		 }
		 $go_split1[$i2]=$key;
		 $go_split2[$i2]=$hash_GO_from_neighbour{$key};
		 $i_neighbour++;
	 }
	 ########### output the potential go terms based on the frequency ########
	 if($i_neighbour == 0)
	 {
		 print "Cannot find any good go potential for $gene_lists[$i], next;\n";
		 next;
	 }
	 print $OUT $gene_lists[$i]."|".$go_split1[0];
	 for($i2=1;($i2<@go_split1) && ($i2<$NUM_prediction);$i2++)
	 {
		 print $OUT "\t".$go_split1[$i2];
	 }
	 print $OUT "\n";

 }



 $OUT->close();



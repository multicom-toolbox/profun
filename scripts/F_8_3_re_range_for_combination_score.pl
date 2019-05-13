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


	  print "This script will re-range the prediction data for the combination score in CAFA 2013, we first convert the score to 2 decimal, and find out the maximum score, then add the maximum score for all scores, such that the final max score is 1. \n";
	  print "perl $0 addr_input addr_output\n";
	  print "For example:\n";
#	  print "perl $0 ../result/processed_map_best_protein_and_including_hitted_protein_by_MHS_and_network_CAFA1 ../test/re_ranged_map_best_protein_including_hitted_protein_MHS_network_CAFA1\n";
#    print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_mining_rules_and_network_CAFA2 ../CAFA_2013_result/F_8_reranged_filtered_valid_prediction_by_MHS_mining_rules_and_network_CAFA2\n";
    print "perl $0 ../CAFA_2013_result/F_8_4_combination_score_CAFA2/Final_combined_file ../CAFA_2013_result/F_8_combination_three_scores_for_CAFA2\n";
          print "perl $0 ../CAFA_2013_result/Test_F_8_2_for_other_species_only/Final_combined_file ../CAFA_2013_result/F_8_combination_three_scores_for_other_species_CAFA2\n";

	  exit(0);
	}

  my($input)=$ARGV[0];
  my($output)=$ARGV[1];

  my($IN,$OUT,$line,$key,$value,$i);
  my(@tem_split,@tem2);
  my(%hash)=();                          # key is target name, value is modelname targetname GO currentscore|mappedscore
  $OUT = new FileHandle ">$output";
  #### 1. First get the maximum score ####
  my($max_score)=-1;
  $IN=new FileHandle "$input";
  while(defined($line=<$IN>))
  {
	  chomp($line);
	  @tem_split=split(/\s+/,$line);

	  if(@tem_split == 4)
	  {#modelname targetname GO score
       $tem_split[3]=sprintf("%.2f",$tem_split[3]);            # convert to 0.01 - 0.99 format
       if($tem_split[3] > $max_score)
       {
       	  $max_score = $tem_split[3];
       }

	  }
	  elsif(@tem_split == 3)
	  {# targetname GO score
		  $tem_split[2]=sprintf("%.2f",$tem_split[2]);            # convert to 0.01 - 0.99 format
       if($tem_split[2] > $max_score)
       {
       	  $max_score = $tem_split[2];
       }	    
	  }
  }
  $IN->close();
  print "We get the maximum score $max_score!\n";
  if($max_score == -1)
  {
  	  print "Check the maximum score, the input file $input has no scores!\n";
  	  exit(0);
  }
  $max_score = 1-$max_score;       # get the difference
  if($max_score < 0)
  {
  	  print "Check the input file, the difference of max score and 1 is negative! $input\n";
  	  exit(0);
  }
  ###  2. Now convert the score to 2 decimal ######
  $IN=new FileHandle "$input";
  while(defined($line=<$IN>))
  {
	  chomp($line);
	  @tem_split=split(/\s+/,$line);

	  if(@tem_split == 4)
	  {#modelname targetname GO score
       $tem_split[3]=sprintf("%.2f",$tem_split[3]);            # convert to 0.01 - 0.99 format
       $tem_split[3]+=$max_score;
       $tem_split[3]=sprintf("%.2f",$tem_split[3]);
       print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\t".$tem_split[3]."\n";
	  }
	  elsif(@tem_split == 3)
	  {# targetname GO score
		  $tem_split[2]=sprintf("%.2f",$tem_split[2]);            # convert to 0.01 - 0.99 format
      $tem_split[2]+=$max_score;
      $tem_split[2]=sprintf("%.2f",$tem_split[2]);
      print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\n";
	  }
  }
  $IN->close();  
  
  $OUT->close();

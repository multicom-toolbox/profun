 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables


   if (@ARGV != 7)
   { # @ARGV used in scalar context = number of args
	  print "This is a latest version at 12/23/2013! We want to keep only two decimals!\n";

	  print "This is more efficient program, should be fast! This program is for CAFA2, combination of three methods! This is for processing large data! \n";

	  print " this program combines the three method's result by giving the weight for each of them! \n";
	  print "perl $0 addr_MIS weight_MIS addr_network weight_network addr_seq_based weight_seq addr_output\n";
	  print "For example:\n";
	  #print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_minig_rules_CAFA2 0.438 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_network_only_CAFA2 0.202 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_sequence_CAFA2 0.251 ../test/old_weight_F_8_comination_three_scores_CAFA2\n";
	  print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_minig_rules_CAFA2 0.876 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_network_only_CAFA2 0.404 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_sequence_CAFA2 0.502 ../test/Two_decimal_New_weight_F_8_comination_three_scores_CAFA2\n";

	  exit(0);
   }

 my($p_MIS)=$ARGV[0];
 my($w_MIS)=$ARGV[1];
 my($p_net)=$ARGV[2];
 my($w_net)=$ARGV[3];
 my($p_seq)=$ARGV[4];
 my($w_seq)=$ARGV[5];
 my($output)=$ARGV[6];

 my(%target2GO)=();           # store the GO function for each target 

 my($IN,$OUT,$line,$i,$j,$new_GO,$key);
 my(@tem_split,@tem222);
 $|=1;
 print "processing $p_MIS ...";
 ######### 1. process MIS ##########
 $IN = new FileHandle "$p_MIS";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {
		 print "Skiping $line!\n";
		 next;
	 }
	 elsif(@tem_split == 3)
	 {# we only have T96060007087    GO:0005737      1.000000.
		 $tem_split[3]=$tem_split[2];
		 $tem_split[2]=$tem_split[1];
		 $tem_split[1]=$tem_split[0];
	 }
	 $tem_split[3]*=$w_MIS;       # update the score using the weight

     $key=$tem_split[1]."_".$tem_split[2];
	 if(not exists $target2GO{$key})
	 {
		 $target2GO{$key} = $tem_split[3];
	 }
	 else
	 {# we already have some GO for this target
		 $target2GO{$key} = 1-(1-$target2GO{$key})*(1-$tem_split[3]);
	 }
 }
 $IN->close();
 print "Done\nprocessing $p_net ...";
 ######### 2. process network based prediction ##########
 $IN = new FileHandle "$p_net";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {
		 print "Skiping $line!\n";
		 next;
	 }
	 elsif(@tem_split == 3)
	 {# we only have T96060007087    GO:0005737      1.000000.
		 $tem_split[3]=$tem_split[2];
		 $tem_split[2]=$tem_split[1];
		 $tem_split[1]=$tem_split[0];
	 }
	 $tem_split[3]*=$w_net;       # update the score using the weight

######## the file is too large, we remove some GO terms which are not very useful ##########
     if($tem_split[3] < 0.005)
	 {
		 next;
	 }
############################################################################################
	 $key=$tem_split[1]."_".$tem_split[2];
	 if(not exists $target2GO{$key})
	 {
		 $target2GO{$key} = $tem_split[3];
	 }
	 else
	 {# we already have some GO for this target
		 $target2GO{$key} = 1-(1-$target2GO{$key})*(1-$tem_split[3]);
	 }
 }
 $IN->close();
 print "Done\n";
=pod 
 my(%hash)=();

 foreach $key (keys %target2GO)
 {
	 $target2GO{$key}=sprintf("%.5f",$target2GO{$key});
	 @tem_split=split(/\_/,$key);
	 if(not exists $hash{$tem_split[0]})
	 {
		 $hash{$tem_split[0]} = $tem_split[1]."_".$target2GO{$key};
	 }
	 else
	 {
		 $hash{$tem_split[0]}.="|".$tem_split[1]."_".$target2GO{$key};
	 }
 }
 %target2GO=();
=cut
 $OUT = new FileHandle ">$output";
 print "processing $p_seq ...";
 ######### 3. process sequence based prediction ##########
 my(@t_add,@tmin);
 
 my($progress)=0;
 my($total_progress)=661141844;

 $IN = new FileHandle "$p_seq";
 while(defined($line=<$IN>))
 {
	 chomp($line);
     $progress++;
	 if($progress%10000 == 0)
	 {
		 print "$progress / $total_progress ... waiting for the program!\n";
	 }

	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {
		 print "Skiping $line!\n";
		 next;
	 }
	 elsif(@tem_split == 3)
	 {# we only have T96060007087    GO:0005737      1.000000.
		 $tem_split[3]=$tem_split[2];
		 $tem_split[2]=$tem_split[1];
		 $tem_split[1]=$tem_split[0];
	 }
	 $tem_split[3]*=$w_seq;       # update the score using the weight
	 $key=$tem_split[1]."_".$tem_split[2];
	 
	 if(not exists $target2GO{$key})
	 {
		 $tem_split[3] = sprintf("%.2f",$target2GO{$key});
		 $target2GO{$key} = $tem_split[3];
	 }
	 else
	 {# we already have some GO for this target
		  $tem_split[3] = 1-(1-$target2GO{$key})*(1-$tem_split[3]);
		  $tem_split[3] = sprintf("%.2f",$tem_split[3]);
		  if($tem_split[3] < 0.01)
		  {
			  $tem_split[3] = 0.01;
		  }
          $target2GO{$key} = $tem_split[3];
	 }
 }
 $IN->close();
 foreach $key (keys %target2GO)
 {
	 $i = sprintf("%.2f",$target2GO{$key});
	 if($i<0.01)
	 {
		 $i=0.01;
	 }
	 $target2GO{$key}=$i;
 }


 ###### output the last target ########
 my($tmp_output)=$output.".tmp";
 $OUT = new FileHandle ">$tmp_output";
 foreach $i (sort {$target2GO{$b} cmp $target2GO{$a}} keys %target2GO)
 {
	 @tem_split=split(/\_/,$i);
	 print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$target2GO{$i}."\n";
 }


 print "Done!\n";
 
 $OUT->close();

 print "Re processing the output file!\n";
 %target2GO=();
 $IN = new FileHandle "$tmp_output";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split!=3)
	 {
		 print "Check $line!!! error in $output, don't touch the output file!\n";
		 exit(0);
	 }
	 if(not exists $target2GO{$tem_split[0]})
	 {
		 $target2GO{$tem_split[0]}=$tem_split[1]."_".$tem_split[2];
	 }
	 else
	 {
		 $target2GO{$tem_split[0]}.="|".$tem_split[1]."_".$tem_split[2];
	 }
 }
 $IN->close();


 my($j);
 $OUT = new FileHandle ">$output";
 foreach $key (keys %target2GO)
 {# $i is the target name
	 @tem_split=split(/\|/,$target2GO{$key});
	 $j=0;
	 for($i=0;$i<@tem_split;$i++)
	 {
		 $j++;
		 if($j>1500)
		 {# we only keep 1500 predictions.
			 last;
		 }
		 @tem222=split(/\_/,$tem_split[$i]);
		 print $OUT "Comined_three_methods\t".$key."\t".$tem222[0]."\t".$tem222[1]."\n";
	 }
	 
 }
 $OUT->close();
 print "Done!\n";


=pod 
 print "processing $p_seq ...";
 ######### 3. process sequence based prediction ##########
 $IN = new FileHandle "$p_seq";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<3)
	 {
		 print "Skiping $line!\n";
		 next;
	 }
	 elsif(@tem_split == 3)
	 {# we only have T96060007087    GO:0005737      1.000000.
		 $tem_split[3]=$tem_split[2];
		 $tem_split[2]=$tem_split[1];
		 $tem_split[1]=$tem_split[0];
	 }
	 $tem_split[3]*=$w_seq;       # update the score using the weight

	 $key=$tem_split[1]."_".$tem_split[2];

######## the file is too large, we remove some GO terms which are not very useful ##########
     if(not exists $target2GO{$key})
	 {# this will definitely not included in the final result, since we keep 0.01 decimals, and this one will be less than 0.01. if already exists in the table, we just calculate that, no waste of memory
		if($tem_split[3] < 0.005)
		{
			 next;
		}
	 }
############################################################################################
	 if(not exists $target2GO{$key})
	 {
		 $target2GO{$key} = $tem_split[3];
	 }
	 else
	 {# we already have some GO for this target
		 $target2GO{$key} = 1-(1-$target2GO{$key})*(1-$tem_split[3]);
	 }
 }
 $IN->close();
 print "Done!\n";
 foreach $key (keys %target2GO)
 {
	 $target2GO{$key}=sprintf("%.5f",$target2GO{$key});
 }
 $OUT = new FileHandle ">$output";
 foreach $i (sort {$target2GO{$b} cmp $target2GO{$a}} keys %target2GO)
 {
	 @tem_split=split(/\_/,$i);
	 print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$target2GO{$i}."\n";
 }

 $OUT->close();
 %target2GO=();
 $IN = new FileHandle "$output";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split!=3)
	 {
		 print "Check $line!!! error in $output, don't touch the output file!\n";
		 exit(0);
	 }
	 if(not exists $target2GO{$tem_split[0]})
	 {
		 $target2GO{$tem_split[0]}=$tem_split[1]."_".$tem_split[2];
	 }
	 else
	 {
		 $target2GO{$tem_split[0]}.="|".$tem_split[1]."_".$tem_split[2];
	 }
 }
 $IN->close();


 my($j);
 $OUT = new FileHandle ">$output";
 foreach $key (keys %target2GO)
 {# $i is the target name
	 @tem_split=split(/\|/,$target2GO{$key});
	 $j=0;
	 for($i=0;$i<@tem_split;$i++)
	 {
		 $j++;
		 if($j>1500)
		 {# we only keep 1500 predictions.
			 last;
		 }
		 @tem222=split(/\_/,$tem_split[$i]);
		 print $OUT "Comined_three_methods\t".$key."\t".$tem222[0]."\t".$tem222[1]."\n";
	 }
	 
 }
 $OUT->close();
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 sub merge_files($$$);     # merge two files
 sub direct_merge_files($$$);

   if (@ARGV != 7)
   { # @ARGV used in scalar context = number of args
	  print "This is a latest version at 12/26/2013! We want to combine three methods' result by first split them!\n";

	  print "This is more efficient program, should be fast! This program is for CAFA2, combination of three methods! This is for processing large data! \n";

	  print " this program combines the three method's result by giving the weight for each of them! \n";
	  print "perl $0 addr_MIS weight_MIS addr_network weight_network addr_seq_based weight_seq addr_output\n";
	  print "For example:\n";
	  #print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_minig_rules_CAFA2 0.438 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_network_only_CAFA2 0.202 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_sequence_CAFA2 0.251 ../test/old_weight_F_8_comination_three_scores_CAFA2\n";
	  print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_minig_rules_CAFA2 0.876 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_network_only_CAFA2 0.404 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_sequence_CAFA2 0.502 ../test/splitting_technique_New_weight_F_8_comination_three_scores_CAFA2\n";
          print "perl $0 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_mining_rules_for_other_species_CAFA2 0.876 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_network_only_for_other_species_CAFA2 0.404 ../CAFA_2013_result/F_7_filtered_valid_prediction_by_sequence_for_other_species_CAFA2 0.502 ../CAFA_2013_result/Test_F_8_2_for_other_species_only\n";


	  exit(0);
   }

 my($p_MIS)=$ARGV[0];
 my($w_MIS)=$ARGV[1];
 my($p_net)=$ARGV[2];
 my($w_net)=$ARGV[3];
 my($p_seq)=$ARGV[4];
 my($w_seq)=$ARGV[5];
 my($dir_output)=$ARGV[6];

 -s $dir_output || system("mkdir $dir_output");
 
 my(%target2GO)=();           # store the GO function for each target 

 my($IN,$OUT,$line,$i,$j,$new_GO,$key);
 my(@tem_split,@tem222);
 $|=1;
 
 my(%all_targets)=();
 my($count_all_targets)=0;
 
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

  if(not exists $all_targets{$tem_split[1]})
  {
  	$all_targets{$tem_split[1]} = 1;
  	$count_all_targets++;
 	}


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
 my(%half_targets)=();
 my($count_half_targets)=0;
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

   if(not exists $half_targets{$tem_split[1]})
   {
   	  if($count_half_targets < $count_all_targets/2)
   	  {
   	  	$half_targets{$tem_split[0]} = 1;
   	  	$count_half_targets++;
   	  }
   	  
   }

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
 
 my($output1)=$dir_output."/"."first_split1";
 my($output2)=$dir_output."/"."first_split2";
 my($OUT1) = new FileHandle ">$output1";
 my($OUT2) = new FileHandle ">$output2";
 
 foreach $key (sort {$target2GO{$b} cmp $target2GO{$a}} keys %target2GO)
 {
 	  @tem_split= split(/\_/,$key);
 	  if(exists $half_targets{$tem_split[0]})
 	  {
 	  	 print $OUT1 $tem_split[0]."\t".$tem_split[1]."\t".$target2GO{$key}."\n";
 	  }
 	  else
 	  {
 	  	 print $OUT2 $tem_split[0]."\t".$tem_split[1]."\t".$target2GO{$key}."\n";
 	  }
 }
 
 
 $OUT1->close();
 $OUT2->close();
 %target2GO=();
 undef %target2GO;
 ######### output the first split end #############
  
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

 $output1=$dir_output."/"."second_split1";
 $output2=$dir_output."/"."second_split2";
 
 $OUT1 = new FileHandle ">$output1";
 $OUT2 = new FileHandle ">$output2";
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
	 
	 if(exists $half_targets{$tem_split[1]})
	 {
	 	  print $OUT1 $tem_split[1]."\t".$tem_split[2]."\t".$tem_split[3]."\n";
	 }
	 else
	 {
	 	  print $OUT2 $tem_split[1]."\t".$tem_split[2]."\t".$tem_split[3]."\n";
	 }
	 
 }
 $OUT1->close();
 $OUT2->close();

 ###### free all memories ######
 %half_targets=();
 undef %half_targets;
 
 ###### combine the targets ########
 print "merging first part ..\n";
 my($output1) = $dir_output."/"."first_split1";
 my($output2) = $dir_output."/"."second_split1";
 my($output3)=$dir_output."/"."Combined_1";
 merge_files($output1,$output2,$output3);
 
print "merging second part!\n";
 $output1 = $dir_output."/"."first_split2";
 $output2 = $dir_output."/"."second_split2";
 my($output3)=$dir_output."/"."Combined_2";
 merge_files($output1,$output2,$output3);
 print "direct merging ..\n"; 
 $output1=$dir_output."/"."Combined_1";
 $output2=$dir_output."/"."Combined_2";
 $output3=$dir_output."/"."Final_combined_file";
 direct_merge_files($output1,$output2,$output3);






  sub direct_merge_files($$$)
  {
  	my($input1,$input2,$output)=@_;
  	my($IN,$OUT,$line);
  	$OUT = new FileHandle ">$output";
  	$IN = new FileHandle "$input1";
  	while(defined($line=<$IN>))
  	{
  		 chomp($line);
  		 print $OUT $line."\n";
  	}
  	$IN->close();
  	$IN = new FileHandle "$input2";
  	while(defined($line=<$IN>))
  	{
  		 chomp($line);
  		 print $OUT $line."\n";
  	}
  	$IN->close();  	
  	$OUT->close();
  }

  sub merge_files($$$)
  {
  	my($input1,$input2,$output)=@_;
  	my($value,$key,$IN,$OUT,$line);
  	my(@tem_split);
  	my(%t2g)=();
  	
  	$IN = new FileHandle "$input1";
  	while(defined($line=<$IN>))
  	{
  		 chomp($line);
  		 @tem_split=split(/\s+/,$line);
  		 if(@tem_split<3)
  		 {
  		 	  print "Warning, check $line\n";
  		 	  next;
  		 }
  	   $key=$tem_split[0]."_".$tem_split[1];
           $tem_split[2] = sprintf("%.4f",$tem_split[2]);
  	   if(not exists $t2g{$key})
  	   {
  	   	  $t2g{$key}=$tem_split[2];
  	   }	
  	   else
  	   {
  	   	  $value = 1-(1-$tem_split[2])*(1-$t2g{$key});
                  $value = sprintf("%.4f",$value);
                  $t2g{$key} = $value;
  	   }
  	}
  	$IN->close();
  	$IN = new FileHandle "$input2";
  	while(defined($line=<$IN>))
  	{
  		 chomp($line);
  		 @tem_split=split(/\s+/,$line);
  		 if(@tem_split<3)
  		 {
  		 	  print "Warning, check $line\n";
  		 	  next;
  		 }
  	   $key=$tem_split[0]."_".$tem_split[1];

           $tem_split[2] = sprintf("%.4f",$tem_split[2]);
  	   if(not exists $t2g{$key})
  	   {
  	          if($tem_split[2] < 0.005)
                  {
                      next;
                  }   
        	  $t2g{$key}=$tem_split[2];
  	   }	
  	   else
  	   {
  	   	  $value = 1-(1-$tem_split[2])*(1-$t2g{$key});
                  $value = sprintf("%.4f",$value);
                  $t2g{$key}=$value;
  	   }
  	}
  	$IN->close();
  	$OUT = new FileHandle ">$output";
    foreach $key (sort {$t2g{$b} cmp $t2g{$a}} keys %t2g)
    {
 	     @tem_split= split(/\_/,$key);  
 	     print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$t2g{$key}."\n";
  	  	
    }
    $OUT->close();  
      %t2g=();  
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
	     if(not exists $t2g{$tem_split[0]})
	     {
		     $t2g{$tem_split[0]}=$tem_split[1]."_".$tem_split[2];
	     }
	     else
	     {
		     $t2g{$tem_split[0]}.="|".$tem_split[1]."_".$tem_split[2];
	     }
    }
    $IN->close();


    my($j);
    my(@tem222);
    $OUT = new FileHandle ">$output";
    foreach $key (keys %t2g)
    {# $i is the target name
	     @tem_split=split(/\|/,$t2g{$key});
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
    %t2g=();
    undef %t2g;
  }
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

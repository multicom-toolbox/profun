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
 use POSIX qw/ceil/;
 sub Mining_rule($);
 if (@ARGV != 6)
 {
         print "This program will evaluate the performance of using different minimum support for apriori to mining the rules based on a training and testing dataset. One GO term in the testing dataset will be randomly removed, and use the rules generated from training dataset to make prediction, and evaluate the performance \n";
		 print "perl $PROGRAM_NAME addr_apriori min_support(0-100) min_confidence(0-100) addr_training_dataset addr_testing_dataset dir_output\n";
		 print "For example:\n";
		 print "perl $PROGRAM_NAME /exports/store1/tool/apriori/src/apriori 0.1 60 ../result/R_2_cross_5_validation_for_GO_terms/training_0 ../result/R_2_cross_5_validation_for_GO_terms/testing_0 ../result/R_3_evaluation_results/s_0.1_c_60_dataset_0_performance\n";
         print "perl $0 /exports/store1/tool/apriori/src/apriori 0.2 60 ../result/R_2_cross_5_validation_for_GO_terms/training_0 ../result/R_2_cross_5_validation_for_GO_terms/testing_0 ../result/R_3_evaluation_results/s_0.2_c_60_dataset_0_performance\n";

	 exit(0);
 }

 my($addr_apri)=$ARGV[0];
 my($min_support)=$ARGV[1];
 my($min_con)=$ARGV[2];
 my($addr_training)=$ARGV[3];
 my($addr_testing)=$ARGV[4];
 my($dir_output)=$ARGV[5];


 -s $dir_output || system("mkdir $dir_output");
 -s $addr_apri || die "Cannot open the C program apriori, $addr_apri\n";
 -s $addr_training || die "Cannot open $addr_training\n";
 -s $addr_testing || die "Cannot open $addr_testing\n";
 if($min_support<0 || $min_support > 100 || $min_con <0 || $min_con>100)
 {
	 print "Check the range of min support and min confidence for using apriori, $min_support, $min_con\n";
	 exit(0);
 }
 $| = 1;
 print "randomly leave one GO term out ...";
 ########################################################################################################################
 #
 #
 #                  1. Randomly take one GO term for each line in the testing dataset, and use that to do evaluation
 #
 #
 ########################################################################################################################
 my($OUT,$OUT2,$line,$IN,$path_in,$value,$path_leaving_one_testing,$path_leaving_one,$ran_one,$key,$i,$j);
 my(@tem_split);
 my($total_for_testing)=0;
 $path_leaving_one_testing = $dir_output."/"."After_leaving_one_GO_term_out_testing_dataset";          # aftering leaving one GO term out for each line, the testing dataset
 $path_leaving_one= $dir_output."/"."The_leaving_out_one_GO_term_for_testing_dataset";        # the leaving one GO terms for the testing dataset
 $OUT = new FileHandle ">$path_leaving_one_testing";
 $OUT2 = new FileHandle ">$path_leaving_one";
 $IN = new FileHandle "$addr_testing";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<1)
	 {
		 next;
	 }
	 elsif(@tem_split == 1)
	 {# only one GO term
		 print $OUT "NULL\n";
		 print $OUT2 $line."\n";
	 }
	 else
	 {
         $total_for_testing++;
	     $ran_one = int(rand(@tem_split));
		 for($i=0;$i<@tem_split;$i++)
		 {
			 if($i == $ran_one)
			 {
				 print $OUT2 $tem_split[$i]."\n";
			 }
			 else
			 {
				 print $OUT $tem_split[$i]."\t";
			 }
		 }
		 print $OUT "\n";
	 }

 }
 $IN->close();
 $OUT->close();
 $OUT2->close();
 print "There are $total_for_testing number of trasactions for testing, leaving one out done!\n";
 print "Generating rules by apriori from training dataset ...";
 ########################################################################################################################
 #
 #
 #                  2. Generate rules from training dataset
 #
 #my($addr_apri)=$ARGV[0];
 #my($min_support)=$ARGV[1];
 #my($min_con)=$ARGV[2];
 #my($addr_training)=$ARGV[3];
 #
 ########################################################################################################################
 my($return_val);
 my($path_rules)=$dir_output."/"."Rules_from_training_dataset";
 $return_val = system("$addr_apri -tr -s$min_support -c$min_con $addr_training $path_rules");
 if($return_val!=0)
 {# we set min support as 0, support means the frequency of that rule exists in the whole database, it's OK if only exists few times, we need high confidence, confidence(x=>y) = support(X and Y)/Support(X). 
	 print "The program :$addr_apri -tr -s$min_support -c$min_con $addr_training $path_rules fails!\n";
 }
 else
 {
	 print "Great, done!\n";
 }
 print "loading rules ...";
 ########################################################################################################################
 #
 #
 #                  3. Loading rules
 #
 #
 ######################################################################################################################## 
   my(%rules)=();                 # store the rules, key is the pre-rule, value is later rule and also include the confidence score, GO:1 => GO:2;100|GO:3;99, so here GO:1 => GO:2 is a rule, with confidence score 100. and GO:1=>GO:3, with conf 99
   my($max_num_com)=0;            # the maximum number of combination elements, for example, if it's 1, then the combination has only one element.
   my(@tem_2,@tem_3);

   $IN=new FileHandle "$path_rules";
   while(defined($line=<$IN>))
   {
	   chomp($line);
	   $line=~s/\s+$//;                # such as :     GO:0004618 <- GO:0006096 GO:0005524 GO:0005737 (0.163133, 82.1343)
       @tem_split=split(/\(/,$line);
	   if(@tem_split<1)
	   {
		   next;
	   }


       @tem_2=split(/\s+/,$tem_split[0]);       #the first is value, rest of them is key

	   if(@tem_2-2 > $max_num_com)
	   {# the max number of combination elements
		   
		   $max_num_com=@tem_2-2;
		   print "Updating at $line, set to $max_num_com\n";
	   }
	   @tem_3=split(/\)/,$tem_split[1]);
	   @tem_3=split(/\s+/,$tem_3[0]);           # this is the confidence score
	   $key=$tem_2[2];
	   for($i=3;$i<@tem_2;$i++)
	   {
		   $key.="_".$tem_2[$i];
	   }
	   $value=$tem_2[0].";".$tem_3[1];
	   if(not exists $rules{$key})
	   {
		   $rules{$key}=$value;
	   }
	   else
	   {
		   $rules{$key}.="|".$value;
	   }
   }


   $IN->close();

   print "The maximum number of combination elements is set to be $max_num_com!  ";

   print "Done!\n";
 
 print "Making prediction by applying the mining rules ...";
 ########################################################################################################################
 #
 #
 #                  4. Use the rules to make prediction
 #
 #
 ########################################################################################################################
 my($path_pre)=$dir_output."/"."Prediction";
 my($prediction,$input_for_pre);
 $OUT = new FileHandle ">$path_pre";
 $IN = new FileHandle "$path_leaving_one_testing";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split<1)
	 {
		 next;
	 }
     if($line eq "NULL")
	 {# there is only one GO term in this line, we skip make prediction
		 print $OUT "SKIP\n";
	 }
	 else
	 {# make prediction 
		 $input_for_pre=$tem_split[0];
		 for($i=1;$i<@tem_split;$i++)
		 {
			 $input_for_pre.="_".$tem_split[$i];
		 }
         $prediction=Mining_rule($input_for_pre);
		 print $OUT $prediction."\n";
	 }
 }
 $IN->close();
 $OUT->close();
 print "Done!\n";
 print "Evaluation ...";
 ########################################################################################################################
 #
 #
 #                  5. Evaluation
 #
 #  1. $path_pre as prediction, SKIP is for only one GO terms, skip it, NULL is for making no prediction, consider it
 #  2. $path_leaving_one, is the one GO term for making prediction
 #
 ########################################################################################################################
 my($total_tp_fp)=0;           # use this and TP to calculate the precision
 my($total_tp_fn)=0;           # use this and TP to calculate the recall
 my($tp)=0;
 my(@pres)=();
 my(@reals)=();
 my($index_for_pre)=0;
 my($index_for_real)=0;
 $IN= new FileHandle "$path_pre";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 $pres[$index_for_pre]=$line;
	 $index_for_pre++;
 }

 $IN->close();
 $IN= new FileHandle "$path_leaving_one";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 $reals[$index_for_real]=$line;
	 $index_for_real++;
 }
 
 $IN->close();
 if($index_for_pre != $index_for_real)
 {
	 die "Check $path_pre and $path_leaving_one, not identical number of elements in these two files, prediction and real functions\n";
 }
 for($i=0;$i<$index_for_pre;$i++)
 {# compare each predictions
	 if($pres[$i] eq "SKIP")
	 {# we don't evaluate this case, since only one GO term in this protein
		 next;
	 }
	 $total_tp_fn+=1;                     # each time we have one GO term as real GO
     
	 if($pres[$i] =~m /$reals[$i]/)
	 {# we find the one we want
		 $tp++;
	 }
	 if($pres[$i] eq "NULL")
	 {# make no prediction for this protein, we don't increase the total_tp_fp
		 next;
	 }
	 @tem_split=split(/\|/,$pres[$i]);
	 $total_tp_fp+=@tem_split;
 }
 my($precision,$recall);
 if($total_tp_fp * $total_tp_fn == 0)
 {
	 print "The total tp+fp is $total_tp_fp, tp+fn is $total_tp_fn, should not be 0!!!\n";
	 exit(0);
 }
 $precision=$tp/$total_tp_fp;
 $recall=$tp/$total_tp_fn;

 my($path_eva)=$dir_output."/"."Performance";
 $OUT=new FileHandle ">$path_eva";
 print $OUT "precision : $precision\n";
 print $OUT "recall : $recall\n";
 $OUT->close();


 sub Mining_rule($)
 {
 	 my($GO2M)=@_;      # the GO function to be mined, such as GO:0005829_GO:0005681_GO:0005685_GO:0005686_GO:0046540_GO:0003723_GO:0000395
	 my($mined_result)="NULL"; 

         my(@tem_split,@tem2,@tem3);
         @tem_split = split(/\_/,$GO2M);
         my($i,$j,$key_rule);
         my($OUT,$IN,$line); 
         my(%hash_rule_GO)=();       # key is the GO term generated from rule, the value is the confidence score
		 my(%hash_tmp_rule)=();

             foreach $key (keys %rules) 
			 {
				 $hash_tmp_rule{$key}=0;
				 for($i=0;$i<@tem_split;$i++)
				 {
					 if($key =~ m/$tem_split[$i]/)
					 {
						 $hash_tmp_rule{$key}++;
					 }
				 }
             }
             foreach $key (keys %hash_tmp_rule)
			 {
				 @tem2=split(/\_/,$key);
				 if(@tem2 == $hash_tmp_rule{$key})
				 {# we mapped a rule
#print "Adding rule : $key=> $rules{$key}\n";
                    @tem2=split(/\|/,$rules{$key});
                    for($i=0;$i<@tem2;$i++)
                    {
                        @tem3=split(/\;/,$tem2[$i]);
                        if(not exists $hash_rule_GO{$tem3[0]})
                        {
#print "adding $tem3[0] with $tem3[1]\n";
                            $hash_rule_GO{$tem3[0]} = $tem3[1]/100;    # get the score
                        }
                        else
                        {
                            $hash_rule_GO{$tem3[0]} = 1- (1-$hash_rule_GO{$tem3[0]})*(1-$tem3[1]/100);
                        }
                    }                     
				 }
             }

         #### write back the rule based prediction ####
         foreach $i (keys %hash_rule_GO)
         {
            if($mined_result eq "NULL")
            {
               $mined_result = $i.";".$hash_rule_GO{$i};
            }
            else
            {
               $mined_result.="|".$i.";".$hash_rule_GO{$i};
            }
         }
         return $mined_result;
 }
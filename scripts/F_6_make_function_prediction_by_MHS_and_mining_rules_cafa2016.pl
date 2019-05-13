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

 sub evalue2score($);
 sub Mining_rule($$$);              # this function will try to use the mining rules to make prediction

  if (@ARGV != 7)
    { # @ARGV used in scalar context = number of args
      print("This program will use and F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl and psiblast result, and also R_1_mining_rules_from_uniprot_by_apriori.pl to get the rules for e-value = 1 proteins. Apply the multiple hits statistics and mining rules method to predict function. \n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir(PSI-BLAST result, *.list file) e_threshold(the e threshold is minmum e-value for mapped protein) addr_real_function_for_each_protein addr_mining_rules_by_apriori addr_java_get_all_combination(one java program to get all combinations of sequence of GO terms) addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print "perl $0 ../CAFA_2012_prediction 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/R_1_mining_rules_from_uniprot/output_s_0.05_c_90 TestSequenceAll ../result/F_6_prediction_by_MHS_and_mining_rules_CAFA1\n";

      print("perl $PROGRAM_NAME ../test/test_CAFA_prediction 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/R_1_mining_rules_from_uniprot/Mining_rules TestSequenceAll ../result/F_6_prediction_by_MHS_and_mining_rules\n\n");

      print "perl $0 ../CAFA_2013_data/other 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/R_1_mining_rules_from_uniprot/Mining_rules TestSequenceAll ../CAFA_2013_result/F_6_prediction_by_MHS_and_mining_rules_for_other_species_CAFA2\n";

	  exit(1) ;
    }
  
   my($dir_psi)=$ARGV[0];
   my($threshold)=$ARGV[1];
   my($addr_real_function)=$ARGV[2];
   my($input_rules)=$ARGV[3];
   my($java_combination)=$ARGV[4];
   my($output)=$ARGV[5];
   my($dir_tmp)=$ARGV[6];

   -s $dir_psi || die "Cannot open the directory $dir_psi!\n";
   -s $addr_real_function || die "Cannot open protein real functions, check the script F_4...\n";

   my($IN,$OUT,$line,$file,$key,$target_name,$n_function,$i,$j,$path_input,$path_out,$type,$value,$evalue);
   my(@files,@tem_split,@tem_2,@tem_3);
   
   if($threshold<0 || $threshold>1)
   {
	   print "Warning, the e-value threshold is set to $threshold!!!!!\n";
   }

   $|=1;
   print "loading all protein/gene real functions ...";
##############################################################################
#
#
#               1. Load all protein / gene's real function
#
#
##############################################################################
   my(%real_function)=();          # key is proteinID/geneID , value is real GO term from uniprot, GO1_GO2...
   $IN=new FileHandle "$addr_real_function";
   while(defined($line = <$IN>))
   {
	   chomp($line);
	   $line=~s/\s+$//;
       @tem_split=split(/\s+/,$line);
	   if(@tem_split<2)
	   {
		   print "Skip $addr_real_function/$line\n";
		   next;
	   }
	   if(not exists $real_function{$tem_split[0]})
	   {
		   $real_function{$tem_split[0]}=$tem_split[1];
	   }
	   else
	   {
		   print "Already exists function for $tem_split[0], check $addr_real_function!\n";
		   next;
	   }
   }

   $IN->close();
   print "Done!\n";

##############################################################################
#
#
#               2. Load the mining rules
#
#
##############################################################################
   my(%rules)=();                 # store the rules, key is the pre-rule, value is later rule and also include the confidence score, GO:1 => GO:2;100|GO:3;99, so here GO:1 => GO:2 is a rule, with confidence score 100. and GO:1=>GO:3, with conf 99
   my($max_num_com)=0;            # the maximum number of combination elements, for example, if it's 1, then the combination has only one element.


my(%hash_test)=();


   print "Loading mining rules ...";
   $IN=new FileHandle "$input_rules";
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
if(not exists $hash_test{@tem_2-2})
{
	$hash_test{@tem_2-2}=1;
}
else
{
	$hash_test{@tem_2-2}++;
}

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


foreach $key (keys %hash_test)
{
	print "$key ... $hash_test{$key}\n";
}

#exit(0);

 
##############################################################################
#
#
#               3. Process each Target ...
#
#
##############################################################################
 $OUT =new FileHandle ">$output";       # output the result 
# my($a,$b);
 print "Processing each target ...\n";
 my($MHS_score);                       # the confidence score for MHS
 my(%hash_all_GO)=();
 my(%hash_rule_GO)=();                 # the GO terms by applying mining rules 
 opendir(DIR,"$dir_psi");
 @files=readdir(DIR);
 foreach $file (@files)
 {
         if($file eq "." || $file eq "..")
         {
                 next;
         }
         @tem_split=split(/\./,$file);
         $target_name=$tem_split[0];
         if(@tem_split<2 || $tem_split[1] ne "list")
         {
                 print "the file name is not *.list format, skip $file !\n";
                 next;
         }
         $path_input=$dir_psi."/".$file;

print "Procssing $path_input !\n";

         %hash_all_GO=();             # store all GO functions and the value is their confidence score
         %hash_rule_GO=();            # store all GO functions mined from apriori method
         ######### 1. Use MHS method to make prediction #########
         $IN=new FileHandle "$path_input";
         while(defined($line=<$IN>))
         {
                 chomp($line);
                 $line=~s/\s+$//;

                 @tem_split=split(/\s+/,$line);
                 if(@tem_split<3)
                 {
                         next;
                 }
                 $evalue=$tem_split[@tem_split-1];
                 if(substr($evalue,0,1) eq "e" || substr($evalue,0,1) eq "E")
                 {
                         $evalue="1".$evalue;
                 }

                 if($evalue>$threshold)
                 {#  the evalue is too large, we don't want this hit
                         next;
                 }
                 $MHS_score=evalue2score($evalue);              # convert evalue to a confidence score
                 @tem_split=split(/\|/,$line);
                 $key=$tem_split[1];             # protein AC
                 if(not exists $real_function{$key})
			     {
					 print "no function for $key, next\n";
					 next;
				 }
				 ######### check whether we want to apply ming rules #####

				 if($evalue==0)
			     {
					 if(not exists $hash_rule_GO{$real_function{$key}})
					 {
						 $hash_rule_GO{$real_function{$key}}="NULL";
					 }
				 }
				 #########################################################
				 @tem_split=split(/\_/,$real_function{$key});
				 for($i=0;$i<@tem_split;$i++)
			     {
					 if(not exists $hash_all_GO{$tem_split[$i]})
					 {
						 $hash_all_GO{$tem_split[$i]}=$MHS_score;
					 }
					 else
					 {
						 $hash_all_GO{$tem_split[$i]}.="|".$MHS_score;
					 }
				 }


         }
         $IN->close();
		 ###### now merge the GO prob in the %hash_all_GO, get the final GO prediction using MHS method ##########
		 foreach $key (keys %hash_all_GO) 
		 {
			 @tem_split=split(/\|/,$hash_all_GO{$key});
			 if(@tem_split<2)
			 {# we only have one probability score
                                 $hash_all_GO{$key}=sprintf("%.6f",$hash_all_GO{$key});
				 next;
			 }
			 #print "original GO prediction : $key($hash_all_GO{$key}), merge to : ";
             $j=1;
			 for($i=0;$i<@tem_split;$i++)
			 {
				 $j=$j*(1-$tem_split[$i]);
			 }
			 $j=1-$j;
                         $j=sprintf("%.6f",$j);
			 $hash_all_GO{$key}=$j;
			 #print "$key($j)\n";
		 }

		 ################################### MHS method finished   ###############################################
         ################################### apply mining rules for the hit protein with evalue = 1 ##############
         my($mined_GO_functions);
		 foreach $key (keys %hash_rule_GO) 
		 {
print "Mining rule $key ..\n";
			 $mined_GO_functions=Mining_rule($key,$max_num_com,$java_combination);               # mining the rules and return mined GO terms,  GO:1;0.9|GO:2;0.8 ..., 0.9 is the confidence score of that rule.
			 if($mined_GO_functions eq "NULL")
			 {# we don't mine any rules
				 next;
			 }
			 @tem_split=split(/\|/,$mined_GO_functions);
			 for($i=0;$i<@tem_split;$i++)
			 {
				 @tem_2=split(/\;/,$tem_split[$i]);
				 if(not exists $hash_all_GO{$tem_2[0]})
				 {
					 $hash_all_GO{$tem_2[0]}=$tem_2[1];
				 }
				 else
				 {
					 $hash_all_GO{$tem_2[0]} = 1- (1-$tem_2[1])*(1-$hash_all_GO{$tem_2[0]});       # merge the GO prediction from the rule and the one from the MHS 
				 }
                                 $hash_all_GO{$tem_2[0]}=sprintf("%.6f",$hash_all_GO{$tem_2[0]});
			 }
		 }

         ################ only keep two decimal ###########
         foreach $key (keys %hash_all_GO)
         {
             $hash_all_GO{$key}=sprintf("%.6f",$hash_all_GO{$key});
         }

         ############## sort the GO score and output the result ############
         foreach $key (sort {$hash_all_GO{$b} cmp $hash_all_GO{$a}} keys %hash_all_GO)
		 {
			 print $OUT $target_name."\t".$key."\t".$hash_all_GO{$key}."\n";
         } 




 }
 $OUT->close();

 sub Mining_rule($$$)
 {
	 my($GO2M,$max_num_com,$java_combination)=@_;      # the GO function to be mined, such as GO:0005829_GO:0005681_GO:0005685_GO:0005686_GO:0046540_GO:0003723_GO:0000395
	 my($mined_result)="NULL"; 
         my($now) = time;
         my($current_id) = $now.$$;
         my($tmp_in)="$dir_tmp/tmp.in.$current_id";
         my($tmp_out)="$dir_tmp/tmp.out.$current_id";
         my(@tem_split,@tem2,@tem3);
         @tem_split = split(/\_/,$GO2M);
         my($i,$j,$key_rule);
         my($OUT,$IN,$line); 
         my(%hash_rule_GO)=();       # key is the GO term generated from rule, the value is the confidence score
		 my(%hash_tmp_rule)=();

      if(@tem_split> 5)
      {
             print "We get more than 5 GO terms for this protein, check $GO2M, the java program maybe very slow, I use another method to do this.\n";
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
print "Adding rule : $key=> $rules{$key}\n";
                    @tem2=split(/\|/,$rules{$key});
                    for($i=0;$i<@tem2;$i++)
                    {
                        @tem3=split(/\;/,$tem2[$i]);
                        if(not exists $hash_rule_GO{$tem3[0]})
                        {
print "adding $tem3[0] with $tem3[1]\n";
                            $hash_rule_GO{$tem3[0]} = $tem3[1]/100;    # get the score
                        }
                        else
                        {
                            $hash_rule_GO{$tem3[0]} = 1- (1-$hash_rule_GO{$tem3[0]})*(1-$tem3[1]/100);
                        }
                    }                     
				 }
             }

     }
     else
	 {# we use java program to process this 

         ### first write all GO terms in a file, and let the java program to get all combinations ####
         $OUT = new FileHandle ">$tmp_in";
         defined($OUT) || die "Cannot open $tmp_in, check whether you can access this path, or you may need to change to another temporary path to store the intermiate result\n";
         for($i=0;$i<@tem_split;$i++)
         {
             print $OUT $tem_split[$i]."\n";
         }
         $OUT->close();
         ####### run the java program to get all combinations #####
         $i=system("java $java_combination $tmp_in $max_num_com > $tmp_out");
         if($i!=0)
         {
             print "Error, java $java_combination $tmp_in $max_num_com fails, check the reason!\n";
             exit(0);
         }

 print "Check rule for $GO2M ...\n";        
         ####### extract the output result, and check whether there is a rule for that, %rules, GO1_GO2 => GO1;100|GO2;99 ... #####
         $IN= new FileHandle "$tmp_out";
         while(defined($line=<$IN>))
         {
            chomp($line);
            @tem_split=split(/\s+/,$line);
print "\n&&&&& processing $line\n";
            if(@tem_split == 1)
            {# only one GO term
                @tem2=split(/\[/,$line);
                @tem3=split(/\]/,$tem2[1]);
print "checking $tem3[0] **** ";
                if(exists $rules{$tem3[0]})
                {
                    @tem2=split(/\|/,$rules{$tem3[0]});
                    for($i=0;$i<@tem2;$i++)
                    {
                        @tem3=split(/\;/,$tem2[$i]);
                        if(not exists $hash_rule_GO{$tem3[0]})
                        {
print "adding $tem3[0] with $tem3[1]\n";
                            $hash_rule_GO{$tem3[0]} = $tem3[1]/100;    # get the score
                        }
                        else
                        {
                            $hash_rule_GO{$tem3[0]} = 1- (1-$hash_rule_GO{$tem3[0]})*(1-$tem3[1]/100);
                        }
                    }
                }
            }
            else
            {# we have more than one GO terms, [GO:1, GO:2, GO:3, ..., GO:n] is the output 
                ### first process the first one ####
                @tem2=split(/\,/,$tem_split[0]);
                @tem3=split(/\[/,$tem2[0]);
				$key_rule=$tem3[1];
				### add the middle GO ###
                for($i=1;$i<@tem_split-1;$i++)
                {
                    @tem2=split(/\,/,$tem_split[$i]);
					$key_rule.="_".$tem2[0];
				}
				### add the last GO ###
                @tem2=split(/\]/,$tem_split[@tem_split-1]);
				$key_rule.="_".$tem2[0];
print "Checking $key_rule ~~~~\n";
				######## check whether we have the rule ##########
                if(exists $rules{$key_rule})
                {
print "adding $key_rule=>$rules{$key_rule}!\n";
                    @tem2=split(/\|/,$rules{$key_rule});
                    for($i=0;$i<@tem2;$i++)
                    {
                        @tem3=split(/\;/,$tem2[$i]);
                        if(not exists $hash_rule_GO{$tem3[0]})
                        {

print "adding $tem3[0]  with $tem3[1]\n";
                            $hash_rule_GO{$tem3[0]} = $tem3[1]/100;    # get the score
                        }
                        else
                        {
                            $hash_rule_GO{$tem3[0]} = 1- (1-$hash_rule_GO{$tem3[0]})*(1-$tem3[1]/100);
                        }
					}
				} # end if               

            } # end else



         }
         $IN->close();

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
         system("rm $tmp_in");
         system("rm $tmp_out");
         return $mined_result;
 }

 sub evalue2score($)
 {######  we use log(e)/100-0.2 to convert the evalue to the range of [0,1], if the result is less than 0, we set to 0, larger than 1, we set to 1.
	  my($evalue_score)=@_;
      my(@tem_split);

      my($value);
      if($evalue_score =~ m/e/)
      {
              @tem_split=split(/e/,$evalue_score);
              if($tem_split[0] eq "")
              {
                      $value=$tem_split[1];
              }
              else
              {
                      $value=$tem_split[1] + (log($tem_split[0])) / (log(10));
              }
      }
      elsif ($evalue_score =~ m/E/)
      {
              @tem_split=split(/E/,$evalue_score);
              if($tem_split[0] eq "")
              {
                      $value=$tem_split[1];
              }
              else
              {
                      $value=$tem_split[1] + (log($tem_split[0])) / (log(10));
              }
      }
      else
      {
              if($evalue_score == 0)
              {
                     return 1;
              }
              $value=(log($evalue_score))/(log(10));
      }
      $value=-1*$value/200-0.01;

	  if($value<0)
	  {
		  $value=0;
	  }
	  if($value>1)
	  {
		  $value=1;
	  }
      return $value;
 }

=pod
We input the *.list file, and also the result of last step, we get the functions for all proteins, we use the multiple hits statistics for making function prediction of each target.
=cut

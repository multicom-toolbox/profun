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

  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args
      print("This program will use and F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl and psiblast result, and also R_1_mining_rules_from_uniprot_by_apriori.pl to get the rules for e-value = 1 proteins. Apply the multiple hits statistics and mining rules method to predict function. \n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir(PSI-BLAST result, *.list file) e_threshold(the e threshold is minmum e-value for mapped protein) addr_real_function_for_each_protein addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME ../test/test_CAFA_prediction 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/F_6_prediction_by_MHS_only\n\n");
          print "perl $0 ../CAFA_2013_data/other 0.01 ../result/F_4_GO_function_for_all_protein_gene ../CAFA_2013_result/F_6_prediction_by_MHS_and_mining_rules_for_other_species_CAFA2\n";


	  exit(1) ;
    }
  
   my($dir_psi)=$ARGV[0];
   my($threshold)=$ARGV[1];
   my($addr_real_function)=$ARGV[2];
   my($output)=$ARGV[3];
   
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
#               2. Process each Target ...
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

         ############## sort the GO score and output the result ############
         foreach $key (sort {$hash_all_GO{$b} cmp $hash_all_GO{$a}} keys %hash_all_GO)
		 {
			 print $OUT $target_name."\t".$key."\t".$hash_all_GO{$key}."\n";
         } 




 }
 $OUT->close();
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

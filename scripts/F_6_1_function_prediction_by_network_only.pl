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


  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
      print("This program will use and F_5_make_function_prediction_for_each_protein_from_network.pl and psiblast result, find the best hit protein, and map that protein to a network, and use the function generanized from the network GO statistics to make function prediction. \n");
      print "revised at 10/18/2013, we multiply the network score with the MHS score, so these two scores are comparable!\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir(PSI-BLAST result, *.list file) addr_mapped_function_for_each_protein addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME ../test/test_CAFA_prediction ../result/F_5_function_prediction_from_network ../test/F_7_function_prediction_from_network\n\n");
      print "perl $0 ../CAFA_2012_prediction ../result/F_5_function_prediction_from_network ../result/F_7_function_prediction_from_network\n";
      print "perl $0 ../CAFA_2013_data/all ../CAFA_2013_result/F_5_function_prediction_from_network ../CAFA_2013_result/F_6_prediction_by_network_only_CAFA2\n";
      print "perl $0 ../CAFA_2013_data/other ../CAFA_2013_result/F_5_function_prediction_from_network ../CAFA_2013_result/F_6_prediction_by_network_only_for_other_species_CAFA2\n";

	  exit(1) ;
    }
  
   my($dir_psi)=$ARGV[0];
   my($addr_network_function)=$ARGV[1];
   my($output)=$ARGV[2];
   
   -s $dir_psi || die "Cannot open the directory $dir_psi!\n";
   -s $addr_network_function || die "Cannot open protein real functions, check the script F_5...\n";

   my($IN,$OUT,$net_function,$select,$line,$file,$target_name,$key,$n_function,$i,$path_input,$type,$value);
   my(@files,@tem_split,@tem_2,@tem_3,@all_species);
   

   $|=1;
   print "loading all protein/gene functions predicted from the network...";
##############################################################################
#
#
#               1. Load all protein / gene's network based predicted function
#
#
##############################################################################
   my(%real_function)=();          # key is proteinID/geneID , value is predicted GO term from network prediction, GO1_score|GO2_score2 .. ; GO1'_score1'|GO2'_score2' ...
   $IN=new FileHandle "$addr_network_function";
   while(defined($line = <$IN>))
   {
	   chomp($line);
	   $line=~s/\s+$//;
       @tem_split=split(/\s+/,$line);
	   if(@tem_split<6)
	   {
		   print "Skip $addr_network_function/$line\n";
		   next;
 	   }
           @all_species = split(/\;/,$tem_split[2]);    # get all species the protein mapped
           if(@all_species == 1)
           {# the protein is mapped to only one network
                $net_function = $tem_split[5];
           }
           else
           {# the protein mapped to multiple networks, we prefer human gene-gene interaction network, and not like the other network.
                $select = 0;
                for($i=0;$i<@all_species;$i++)
                {
                       if($all_species[$i] eq "3_contact_NORMAL_B_gene")
                       { 
                           $select = $i;
                           last;
                       }
                       elsif(($all_species[$i] eq "other") && ($select == $i))
                       {# we don't prefer other network
                           $select++;
                       } 
                } 
                if($select >= @all_species)
                { 
                       print "Should not happen, check here, select is $select, and all speicies is ".@all_species."\n";
                       $select = 0;
                }
                @tem_2=split(/\;/,$tem_split[5]);
                $net_function=$tem_2[$select];    # we assign the selected network predicted function
	    }
           if(not exists $real_function{$tem_split[0]})
           {
	       $real_function{$tem_split[0]}=$net_function;;
           }
           else
           {
	       print "Already exists function for $tem_split[0], check $addr_network_function!\n";
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
 my($MHS_score);               # we use the MHS score as the confidence score of that network
 my(%hash_all_GO)=();
 my($evalue);
 print "Processing each target ...\n";
 $OUT = new FileHandle ">$output";          # output the result 
  
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
print "#########processing $path_input ...\n";
         $net_function = "NULL";
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
                 $MHS_score=evalue2score($evalue);              # convert evalueto a confidence score
                 @tem_split=split(/\|/,$line);        # now we need to find the protein AC, tem_split[1]
                 if(not exists $real_function{$tem_split[1]})
                 {# we don't find the best hit protein's function by the network method, we continue to search the second match protein.
#print "We don't get any network based function for $tem_split[1] with evalue $MHS_score\n";
                         next;
                 }
                 elsif($real_function{$tem_split[1]} eq "NULL")
                 {# there is no function for the best hit protein by network method, we continue to find the second match one
#print "There is no function of network based for $tem_split[1] with evalue $MHS_score\n";
                        next;
                 }
                 else
                 {# we find some function by network to make prediction
#print "We get the function for $tem_split[1]";
                        $net_function=$real_function{$tem_split[1]};
#print "See: $net_function\n";
                        last;        # finish searching
                 }
         }
         $IN->close();
         ############ output the predictions #########
         if($net_function eq "NULL")
         {
                print "No function prediction from network can be made for this target : $file, at $path_input, check the reason!\n";
                next;
         }
         %hash_all_GO=();
         @tem_split=split(/\|/,$net_function);
         for($i=0;$i<@tem_split;$i++)
         {
             @tem_2=split(/\_/,$tem_split[$i]);
             
             if(not exists $hash_all_GO{$tem_2[0]})
             {
                 $tem_2[1] = $tem_2[1] * $MHS_score;            # times the MHS score, so now they are comparable
                 $tem_2[1] = sprintf("%.6f",$tem_2[1]);
                 $hash_all_GO{$tem_2[0]} = $tem_2[1];
             }
             else
             {
                print "Should not happen, $tem_2[0] exists more than one time in the function network based prediction, check $net_function also!\n";
                
             }
         } 
         foreach $key (sort {$hash_all_GO{$b} cmp $hash_all_GO{$a}} keys %hash_all_GO)
         {
             print $OUT $target_name."\t".$key."\t".$hash_all_GO{$key}."\n";
         }

 }
 $OUT->close();
 sub evalue2score($)
 {######  we use log(e)/100-0.2 to convert the evalue to the range of [0,1], ifthe result is less than 0, we set to 0, larger than 1, we set to 1.
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
      





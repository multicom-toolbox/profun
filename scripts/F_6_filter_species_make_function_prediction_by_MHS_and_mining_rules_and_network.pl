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

  if (@ARGV != 8)
    { # @ARGV used in scalar context = number of args
      print("This program will first get the species information, and try to filter the BLAST output by filtering the species information,  use and F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl and psiblast result, and also R_1_mining_rules_from_uniprot_by_apriori.pl to get the rules for e-value = 1 proteins, in addition, we use F_5_make_function_prediction_for_each_protein_from_network.pl and psiblast result, find the best hit protein, and map that protein to a network, and use the function generanized from the network GO statistics to make function prediction. Apply the multiple hits statistics and mining rules method, and networks to predict function. \n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/29/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_uniprot_dat dir(PSI-BLAST result, *.list file) e_threshold(the e threshold is minmum e-value for mapped protein) addr_real_function_for_each_protein addr_mining_rules_by_apriori addr_java_get_all_combination(one java program to get all combinations of sequence of GO terms) addr_mapped_function_for_each_protein addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print "perl $0 /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../CAFA_2012_prediction 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/R_1_mining_rules_from_uniprot/output_s_0.05_c_90 TestSequenceAll ../result/F_5_function_prediction_from_network ../result/F_6_prediction_by_filtering_species_and_MHS_and_mining_rules_and_network_CAFA1\n";

      print("perl $PROGRAM_NAME /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../test/test_CAFA_prediction 0.01 ../result/F_4_GO_function_for_all_protein_gene ../result/R_1_mining_rules_from_uniprot/Mining_rules TestSequenceAll ../result/F_5_function_prediction_from_network ../result/F_6_prediction_by_filtering_and_MHS_and_mining_and_network_rules\n\n");

	  exit(1) ;
    }

   my($addr_uniprot)=$ARGV[0];
   my($dir_psi)=$ARGV[1];
   my($threshold)=$ARGV[2];
   my($addr_real_function)=$ARGV[3];
   my($input_rules)=$ARGV[4];
   my($java_combination)=$ARGV[5];
   my($addr_network_function)=$ARGV[6];
   my($output)=$ARGV[7];
   
   -s $dir_psi || die "Cannot open the directory $dir_psi!\n";
   -s $addr_real_function || die "Cannot open protein real functions, check the script F_4...\n";
   -s $addr_network_function || die "Cannot open protein real functions, check the script F_5...\n";

   my($IN,$OUT,$line,$file,$key,$target_name,$n_function,$net_function,$select,$i,$j,$path_input,$path_out,$type,$value,$evalue);
   my(@files,@tem_split,@tem_2,@tem_3,@all_species);
   
   if($threshold<0 || $threshold>1)
   {
	   print "Warning, the e-value threshold is set to $threshold!!!!!\n";
   }

   $|=1;
print "Loading uniprot database and get the species for each protein ... ";
#sleep(1);
#############################################################################################
#
#                First. Load uniprot database
#
#############################################################################################
my(@all_ids);
$all_ids[0] = "186497";
$all_ids[1] = "243232";
$all_ids[2] = "273057";
$all_ids[3] = "309800";
$all_ids[4] = "436308";
$all_ids[5] = "453591";
$all_ids[6] = "478009";
$all_ids[7] = "160488";
$all_ids[8] = "170187";
$all_ids[9] = "208964";

$all_ids[10] = "223283";
$all_ids[11] = "224308";
$all_ids[12] = "243273";
$all_ids[13] = "321314";
$all_ids[14] = "83333";
$all_ids[15] = "85962";
$all_ids[16] = "99287";
$all_ids[17] = "10090";
$all_ids[18] = "10116";
$all_ids[19] = "284812";
$all_ids[20] = "3702";
$all_ids[21] = "44689";
$all_ids[22] = "559292";
$all_ids[23] = "7227";
$all_ids[24] = "7955";
$all_ids[25] = "8355";
$all_ids[26] = "9606";
my(%all_possible_spe)=();     # get all possible species
for($i=0;$i<@all_ids;$i++)
{
    $all_possible_spe{$all_ids[$i]}=1;
}
my($spe_id);

my(%ID2spe)=();              # the protein ID mapped to species


my($index_ID,$index_AC,$key_spe,$i1,$i2);
my(@v_AC,@species_ID,@tem_for_spe);
my(%tem_spe);

@species_ID=();              # store all species ID for this
@v_AC=();
$index_ID=0;
$index_AC=0;

$IN=new FileHandle "$addr_uniprot";
while(defined($line=<$IN>))
{
    chomp($line);
    $line=~s/\s+$//;
    @tem_split=split(/\s+/,$line);
    if(substr($line,0,2) eq "//")
    {# new line comes
        if(@species_ID == 0)
        {
            print "Warning, check this @v_AC, cannot find species for this protein record, check!\n";
            next;
        }
        
        for($i=0;$i<@v_AC;$i++)
        {# build relationship of AC and species
            if(not exists $ID2spe{$v_AC[$i]})
            {
                $key_spe=$species_ID[0];
                for($i1=1;$i1<@species_ID;$i1++)
                {
                    $key_spe.="|".$species_ID[$i1];
                }
                $ID2spe{$v_AC[$i]}=$key_spe;
            }
            else
            {
                print "!!!!!! this protein AC $v_AC[$i] mapped to more than one species ID !!!!\n";
                @tem_for_spe=split(/\|/,$ID2spe{$v_AC[$i]});
                %tem_spe=();
                for($i1=0;$i1<@tem_for_spe;$i1++)
                {
                    if(not exists $tem_spe{$tem_for_spe[$i1]})
                    {
						$tem_spe{$tem_for_spe[$i1]}=1;
                    }
                }
                for($i1=0;$i1<@species_ID;$i1++)
                {
                    if(not exists $tem_spe{$species_ID[$i1]})
                    {
						$tem_spe{$species_ID[$i1]}=1;
                    }
                }
                $j=0;
                foreach $value (keys %tem_spe)
                {
                    if($j==0)
                    {
                        $ID2spe{$v_AC[$i]}=$value;
                        $j=1;
                    }
                    else
                    {
                        $ID2spe{$v_AC[$i]}.="|".$value;
                    }
                }
            }# end else
            
            
        }
        
        
        @species_ID=();              # store all species ID for this
        @v_AC=();
        $index_ID=0;
        $index_AC=0;
        
        next;
    }
    if($tem_split[0] eq "AC")
    {
        for($i=1;$i<@tem_split;$i++)
        {
            @tem_2=split(/\;/,$tem_split[$i]);
            $v_AC[$index_AC]=$tem_2[0];
            $index_AC++;
        }
    }
    elsif(($tem_split[0] eq "OX"))           # here we only consider OX
    {# get the species, such as line is : OX   NCBI_TaxID=654924;
        
        @tem_for_spe=split(/\;/,$tem_split[1]);
        @tem_for_spe=split(/\=/,$tem_for_spe[0]);
        
        $species_ID[$index_ID] = $tem_for_spe[1];
        $index_ID++;
    }
    
    
}

$IN->close();
print "Done!\n";


   print "loading all protein/gene functions predicted from the network...";
##############################################################################
#
#
#               0. Load all protein / gene's network based predicted function
#
#
##############################################################################
   my(%network_function)=();          # key is proteinID/geneID , value is predicted GO term from network prediction, GO1_score|GO2_score2 .. ; GO1'_score1'|GO2'_score2' ...
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
           if(not exists $network_function{$tem_split[0]})
           {
	            $network_function{$tem_split[0]}=$net_function;;
           }
           else
           {
	       print "Already exists function for $tem_split[0], check $addr_network_function!\n";
	       next;
           }
           
   }

   $IN->close();
   print "Done!\n";


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
 my($MHS_score,$MHS_score_network);                       # the confidence score for MHS, and the one for the network used protein
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
        ####### get the taxid of this file #######
        $i=length($tem_split[0])-8;
        $spe_id=substr($tem_split[0],1,$i);         # get the species name
        if(not exists $all_possible_spe{$spe_id})
        {
            print "Fails, check the file $file, the species ID is $spe_id, not possible!\n";
            exit(0);
        }
        ##########################################
     
     
         $net_function = "NULL";

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
             ##################################################
             ###### check whether this protein is in the same species with the query protein. If it is , we it's in the same species, we give the weight 0.8, and if it's not in the same species, we times 0.2.
                 if(not exists $ID2spe{$key} || $ID2spe{$key} ne $spe_id)
                 {
                     $MHS_score=$MHS_score*0.2;
                 }
                 else
                 {
                     $MHS_score=$MHS_score*0.8;
                 }
             ##################################################
             
             
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
                 ##########################################################
				 ##!!!  now get the function predicted from the network !!!###
				 @tem_split=split(/\|/,$line);
				 if($net_function ne "NULL")
			     {# we already find the best match protein's network based funciton prediction, we don't need to continue searching for other's.
					 next;
				 }
                 if(not exists $network_function{$tem_split[1]})
                 {# we don't find the best hit protein's function by the network method, we continue to search the second match protein.
#print "We don't get any network based function for $tem_split[1] with evalue $MHS_score\n";
                         next;
                 }
                 elsif($network_function{$tem_split[1]} eq "NULL")
                 {# there is no function for the best hit protein by network method, we continue to find the second match one
#print "There is no function of network based for $tem_split[1] with evalue $MHS_score\n";
                        next;
                 }
                 else
                 {# we find some function by network to make prediction
#print "We get the function for $tem_split[1]";
                        $net_function=$network_function{$tem_split[1]};
						$MHS_score_network=$MHS_score;
#print "See: $net_function\n";
                 }

         }
         $IN->close();
		 ###### now merge the GO prob in the %hash_all_GO, get the final GO prediction using MHS method ##########
		 foreach $key (keys %hash_all_GO) 
		 {
			 @tem_split=split(/\|/,$hash_all_GO{$key});
			 if(@tem_split<2)
			 {# we only have one probability score
                                 $hash_all_GO{$key}=sprintf("%.2f",$hash_all_GO{$key});
				 next;
			 }
			 #print "original GO prediction : $key($hash_all_GO{$key}), merge to : ";
             $j=1;
			 for($i=0;$i<@tem_split;$i++)
			 {
				 $j=$j*(1-$tem_split[$i]);
			 }
			 $j=1-$j;
                         $j=sprintf("%.2f",$j);
			 $hash_all_GO{$key}=$j;
			 #print "$key($j)\n";
		 }

		 ################################### MHS method finished   ###############################################
         ################################### apply mining rules for the hit protein with evalue = 1 ##############
         my($mined_GO_functions);
		 foreach $key (keys %hash_rule_GO) 
		 {
#print "Mining rule $key ..\n";
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
                 $hash_all_GO{$tem_2[0]}=sprintf("%.2f",$hash_all_GO{$tem_2[0]});
			 }
		 }

         ################################ apply network method and add the function predicted from network #########
		 if($net_function ne "NULL")
	     {# we need the protein has network based function prediction
			 
           @tem_split=split(/\|/,$net_function);
           for($i=0;$i<@tem_split;$i++)
           {
             @tem_2=split(/\_/,$tem_split[$i]);
             $tem_2[1] = $tem_2[1] * (1-$MHS_score_network);            # times the 1- MHS score, so now they are comparable
             $tem_2[1] = sprintf("%.2f",$tem_2[1]);             
             if(not exists $hash_all_GO{$tem_2[0]})
             {
                 $hash_all_GO{$tem_2[0]} = $tem_2[1];
             }
             else
             {
                #print "Should not happen, $tem_2[0] exists more than one time in the function network based prediction, check $net_function also!\n";
				$value=1- (1-$tem_2[1])*(1-$hash_all_GO{$tem_2[0]});       # merge the GO prediction from rule, MHS and the one from the network
				$value=sprintf("%.2f",$value);
                $hash_all_GO{$tem_2[0]} = $value;
             }
           } 
		 }

         ############## sort the GO score and output the result ############
         foreach $key (sort {$hash_all_GO{$b} cmp $hash_all_GO{$a}} keys %hash_all_GO)
		 {
			 if($hash_all_GO{$key} == 0)
			 {
				 next;
			 }
			 print $OUT $target_name."\t".$key."\t".$hash_all_GO{$key}."\n";
         } 




 }
 $OUT->close();

 sub Mining_rule($$$)
 {
	 my($GO2M,$max_num_com,$java_combination)=@_;      # the GO function to be mined, such as GO:0005829_GO:0005681_GO:0005685_GO:0005686_GO:0046540_GO:0003723_GO:0000395
	 my($mined_result)="NULL"; 
         my($tmp_in)="/tmp/tmp.in";
         my($tmp_out)="/tmp/tmp.out";
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

 #print "Check rule for $GO2M ...\n";        
         ####### extract the output result, and check whether there is a rule for that, %rules, GO1_GO2 => GO1;100|GO2;99 ... #####
         $IN= new FileHandle "$tmp_out";
         while(defined($line=<$IN>))
         {
            chomp($line);
            @tem_split=split(/\s+/,$line);
#print "\n&&&&& processing $line\n";
            if(@tem_split == 1)
            {# only one GO term
                @tem2=split(/\[/,$line);
                @tem3=split(/\]/,$tem2[1]);
#print "checking $tem3[0] **** ";
                if(exists $rules{$tem3[0]})
                {
                    @tem2=split(/\|/,$rules{$tem3[0]});
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
#print "Checking $key_rule ~~~~\n";
				######## check whether we have the rule ##########
                if(exists $rules{$key_rule})
                {
#print "adding $key_rule=>$rules{$key_rule}!\n";
                    @tem2=split(/\|/,$rules{$key_rule});
                    for($i=0;$i<@tem2;$i++)
                    {
                        @tem3=split(/\;/,$tem2[$i]);
                        if(not exists $hash_rule_GO{$tem3[0]})
                        {

#print "adding $tem3[0]  with $tem3[1]\n";
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

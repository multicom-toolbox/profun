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
 sub make_prediction_MHS($$);
 sub convert2prob($);        # convert the frequency to prob

  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args
      print("This program will use F_3_get_one_step_neighbours_for_all_mapped_ID.pl and F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl result, use the interaction network GO statistics to make prediction for each protein. \n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/6/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME protein_list(proteinID evalue_target_name species mapped_protein/gene neighbours_of_mapped_protein/gene one_step_neighbours) protein_gene_real_function dir_network_GO_statistics addr_output\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME ../result/F_3_protein_lists_with_network_and_one_step_neighbour ../result/F_4_GO_function_for_all_protein_gene ../result/1_go_frequency ../result/F_5_function_prediction_from_network\n\n");

	  exit(1) ;
    }
  
   my($addr_protein_list)=$ARGV[0];
   my($addr_real_function)=$ARGV[1];
   my($dir_network_stat)=$ARGV[2];
   my($output)=$ARGV[3];
   
   -s $dir_network_stat || die "Cannot open the directory $dir_network_stat, Check script A_3_cal_go_pairs_frequency.pl, B_3_cal_go_pairs_frequency.pl!\n";
   -s $addr_protein_list || die "Cannot open protein lists file, check F_3_get_one_step_neighbours_for_all_mapped_ID.pl\n";

   my($IN,$OUT,$line,$file,$key,$n_function,$i,$path_input,$type);
   my(@files,@tem_split,@tem_2);
   
   $|=1;
   print "loading all protein/gene real functions ...";
##############################################################################
#
#
#               1. Load all protein / gene's real function
#
#
##############################################################################
   my(%real_function)=();          # key is proteinID/geneID , value is real GO term from uniprot
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
#               2. Load GO statistics 
#
#
##############################################################################
   print "Loading GO statistics ...";
   my(%hash_GO_stat)=();                      # key is type_GO, value is the statistics for the GO
   opendir(DIR,"$dir_network_stat");
   @files=readdir(DIR);
   foreach $file (@files)
   {
	   if($file eq "." || $file eq "..")
	   {
		   next;
	   }
	   $path_input=$dir_network_stat."/".$file;
	   @tem_split=split(/\./,$file);
	   if(@tem_split<2)
	   {
		   print "We skip the file for GO statistics $file, since the name is not *.*\n";
		   next;
	   }
	   $type=$tem_split[0];
	   $IN=new FileHandle "$path_input";
	   while(defined($line=<$IN>))
	   {
		   chomp($line);
		   $line=~s/\s+$//;
		   @tem_split=split(/\s+/,$line);
		   if(@tem_split<2)
		   {
			   print "$path_input, We skip $line!\n";
			   next;
		   }
		   $key=$type."|".$tem_split[0];
		   if(not exists $hash_GO_stat{$key})
		   {
			   $hash_GO_stat{$key}=convert2prob($tem_split[1]);
		   }
	   }
	   $IN->close();
   }
   print "Done!\n";
##############################################################################
#
#
#               3. Load protein/gene lists and make prediction
#
#
##############################################################################  
   print "Loading protein/gene lists and make function prediction from the network GO statistics ...\n";
   my(@species)=();
   my(@neighbours)=();
   my(@prediction)=();
   
   $OUT = new FileHandle ">$output";
   $IN = new FileHandle "$addr_protein_list";
   while(defined($line=<$IN>))
   {
	   chomp($line);
	   @tem_split=split(/\s+/,$line);
	   if(@tem_split<5)
	   {
		   print "We need the protein list file with five columns,proteinID evalue_target_name species mapped_protein/gene neighbours_of_mapped_protein/gene one_step_neighbours, check $line in $addr_protein_list\n";
		   next;
	   }
	   if($tem_split[2] eq "NULL" || $tem_split[3] eq "NULL" || $tem_split[4] eq "NULL")
	   {
		   print $OUT $line."\t"."NULL"."\n";
		   next;
	   }
	   @species=split(/\;/,$tem_split[2]);
	   @neighbours=split(/\|/,$tem_split[4]);
	   
	   @prediction=();
       $prediction[0]=make_prediction_MHS($species[0],$neighbours[0]);
	   print $OUT $line."\t".$prediction[0];
	   
	   for($i=1;$i<@species;$i++)
	   {
		   $prediction[$i]=make_prediction_MHS($species[$i],$neighbours[$i]);
		   print $OUT ";".$prediction[$i];
	   }
	   print $OUT "\n";
       

   }
   $IN->close();
   $OUT->close();

  sub make_prediction_MHS($$)
  {
	  my($spe,$nei)=@_;

	  my(%hash_all_GO)=();  # key is GO term,  value is prob1|prob2 ...
      my($my_pre);
	  my(@all_n)=split(/\_/,$nei);
	  my(@GO_fun,@temtem,@temtem2);
	  my($j,$i,$k,$key);
	  for($i=0;$i<@all_n;$i++)
	  {
		  if(not exists $real_function{$all_n[$i]})
		  {
			  print "no function for $all_n[$i]\n";
			  next;
		  }
          $key=$real_function{$all_n[$i]};
		  @GO_fun=split(/\_/,$key); 
		  for($j=0;$j<@GO_fun;$j++)
		  {
			  $key=$spe."|".$GO_fun[$j];
              if(not exists $hash_GO_stat{$key})
			  {
				  print "No GO statistics for $key!\n";
				  next;
			  }
			  @temtem=split(/\|/,$hash_GO_stat{$key});
			  for($k=0;$k<@temtem;$k++)
			  {
				  @temtem2=split(/\_/,$temtem[$k]);
				  if(not exists $hash_all_GO{$temtem2[0]})
				  {
					  $hash_all_GO{$temtem2[0]}=$temtem2[1];
				  }
				  else
				  {
					  $hash_all_GO{$temtem2[0]}.="|".$temtem2[1];
				  }
			  }

		  }
	  }

	  $my_pre="NULL";
	  foreach $key (keys %hash_all_GO) 
	  {
		  @temtem=split(/\|/,$hash_all_GO{$key});
		  $j=1;
		  for($i=0;$i<@temtem;$i++)
		  {
			  $j=$j*(1-$temtem[$i]);
		  }
		  $j=1-$j;
		  if($my_pre eq "NULL")
		  {
			  $my_pre=$key."_".$j;
		  }
		  else
		  {
			  $my_pre.="|".$key."_".$j;
		  }
		  
	  }
	  return $my_pre;


  }

  sub convert2prob($)
  {
	  my($fre)=@_;
	  my(@all_GO)=();
	  @all_GO=split(/\|/,$fre);
	  my($i);
	  my(@tem_split);
	  my($total)=0;
	  my(@GO_list)=();
	  my(@GO_pro)=();
	  for($i=0;$i<@all_GO;$i++)
	  {
		  @tem_split=split(/\_/,$all_GO[$i]);
		  $total+=$tem_split[1];
		  $GO_list[$i]=$tem_split[0];
		  $GO_pro[$i]=$tem_split[1];
	  }
	  if($total==0)
	  {
		  print "Should not happen, check the $fre!\n";
		  exit(0);
	  }
	  my($pro_list);
	  for($i=0;$i<@GO_list;$i++)
	  {
		  $GO_pro[$i]/=$total;
	  }

	  $pro_list=$GO_list[0]."_".$GO_pro[0];
	  for($i=1;$i<@GO_list;$i++)
	  {
		  $pro_list.="|".$GO_list[$i]."_".$GO_pro[$i];
	  }
	  return $pro_list;
  }


#We use the neighbour's function and make predictions

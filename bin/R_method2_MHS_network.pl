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



  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
      print "This program takes the BLAST *.list file folder as input, and make function prediction using MHS method!\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 8/11/2014\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_script(the path of script folder) dir_CAFA_data(folder to store the data for using) dir_BLAST_result e_value_threshold dir_output!\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME \n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";

 my($dir_base) = $ARGV[0];
 my($dir_base_data) = $ARGV[1];
 my($dir_input_BLAST)=$ARGV[2];
 my($threshold)=$ARGV[3];
 my($dir_output)=$ARGV[4];

 -s $dir_input || die "No input file $dir_input!\n";
 -s $dir_output || system("mkdir $dir_output");
  
 my($return_val,$key,$value,$IN,$OUT,$path_input,$target_name,$path_output,$line,$file,$i,$j,$evalue);
 my(@tem_split,@tem_2,@files);
 

 ########### 1. call F_1_get_all_proteins_from_blast_result.pl #############
 my($F_1) = $dir_base."/"."F_1_get_all_proteins_from_blast_result.pl";
 -s $F_1 || die "Cannot find script $F_1, check the path setting\n";
 my($F_1_out) = $dir_output."/"."F_1_out";
 $return_val = system("perl $F_1 $dir_input_BLAST 0.01 $F_1_out");
 if($return_val)
 {
    print "perl $F_1 $dir_input_BLAST 0.01 $F_1_out fails!\n";
    exit(0);
 }
 ########## 2. call F_2_assign_network_type_for_each_protein.pl ############
 my($F_2) = $dir_base."/"."F_2_assign_network_type_for_each_protein.pl";
 -s $F_2 || die "Cannot find script $F_2, check the path setting\n"; 
 my($uniprot) = $dir_base_data."/"."uniprot_sprot.dat";
 my($F_2_net) = $dir_base_data."/"."data_for_F_2_all_networks";
 my($F_2_out) = $dir_output."/"."F_2_out";
 $return_val = system("perl $F_2 $uniprot $F_2_net $F_1_out $F_2_out");
 if($return_val)
 {
    print "perl $F_2 $uniprot $F_2_net $F_1_out $F_2_out fails!\n";
    exit(0);
 } 

 ######### 3. call F_3_get_one_step_neighbours_for_all_mapped_ID.pl#########
 my($F_3) = $dir_base."/"."F_3_get_one_step_neighbours_for_all_mapped_ID.pl";
 -s $F_3 || die "Cannot find script $F_3, check the path setting\n"; 
 my($one_step) = $dir_base_data."/"."3_one_step_protein_neighbours";
 my($F_3_out) = $dir_output."/"."F_3_out";
 $return_val = system("perl $F_3 $one_step $F_2_out $F_3_out");
 if($return_val)
 {
    print "perl $F_3 $one_step $F_2_out $F_3_out fails!\n";
    exit(0);
 }
 
 ######### 4. F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl ######
 my($F_4) = $dir_base."/"."F_4_get_GO_function_for_all_neighbours_and_mapped_ID.pl";
 -s $F_4 || die "Cannot find script $F_4, check the path setting\n";
 my($F_4_out) = $dir_output."/"."F_4_out";
 $return_val = system("perl $F_4 $uniprot $F_3_out $F_4_out");
 if($return_val)
 {
    print "perl $F_4 $uniprot $F_3_out $F_4_out fails!\n";
    exit(0);
 } 

 ######## call F_5_make_function_prediction_for_each_protein_from_network.pl ##
 my($F_5) = $dir_base."/"."F_5_make_function_prediction_for_each_protein_from_network.pl";
 -s $F_5 || die "Cannot find script $F_5, check the path setting\n";
 my($F_5_out) = $dir_output."/"."F_5_out";
 my($go_fre) = $dir_base_data."/"."1_go_frequency";
 $return_val = system("perl $F_5 $F_3_out $F_4_out $go_fre $F_5_out");
 if($return_val)
 {
    print "perl $F_5 $F_3_out $F_4_out $go_fre $F_5_out fails!\n";
    exit(0);
 }

 



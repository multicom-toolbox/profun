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
 use Cwd;
 use Cwd 'abs_path';


  if (@ARGV != 5)
    { # @ARGV used in scalar context = number of args
      print "This program takes the BLAST *.list file folder as input, and make function prediction using three methods:\n";
      print "1. Combination of MHS, network, and sequence based function prediction\n";
      print "2. MHS and applying mining rules for function prediction\n";
      print "3. MHS, network based function prediction\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 8/11/2014\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_script(the path of script folder) dir_CAFA_data(folder to store the data for using) addr_fasta_seq dir_BLAST_result dir_output!\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
#      print("perl $PROGRAM_NAME ../script/ ../CAFA_data/ ../test/test_fasta ../test/BLAST_all_data ../test/R_all_output\n\n");
          print "perl $0 /var/www/html/profunc/CAFA_website_2014/script /var/www/html/profunc/CAFA_website_2014/CAFA_data ../test/test_fasta ../test/BLAST_all_data ../test/R_all_output_test\n";
	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";

 my($dir_base) = abs_path($ARGV[0]);
 my($dir_base_data) = abs_path($ARGV[1]);
 my($addr_fasta) = abs_path($ARGV[2]);
 my($dir_input_BLAST)= abs_path($ARGV[3]);
 my($dir_output)= abs_path($ARGV[4]);

 -s $dir_output || system("mkdir $dir_output");
  
 my($return_val,$key,$value,$IN,$OUT,$path_input,$target_name,$path_output,$line,$file,$i,$j,$evalue);
 my(@tem_split,@tem_2,@files);

 ###### set the environment ########
 $ENV{'CLASSPATH'}.=':$dir_base';
 ####### go to the dir_base folder ###########
 chdir($dir_base);

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

 ###### call F_6_make_function_prediction_by_MHS_and_mining_rules.pl      ####
 my($F_6_MHS) = $dir_base."/"."F_6_make_function_prediction_by_MHS_and_mining_rules.pl";
 -s $F_6_MHS || die "Cannot find script $F_6_MHS, check the path setting\n";
 my($mining_rule) = $dir_base_data."/"."output_s_0.05_c_90";
 my($testsequenceall) = "TestSequenceAll";
 my($F_6_MHS_out) = $dir_output."/"."F_6_MHS_out";
 $return_val = system("perl $F_6_MHS $dir_input_BLAST 0.01 $F_4_out $mining_rule $testsequenceall $F_6_MHS_out");
 if($return_val)
 {
    print "perl $F_6_MHS $dir_input_BLAST 0.01 $F_4_out $mining_rule $testsequenceall $F_6_MHS_out fails!\n";
    exit(0);
 }
 
 ##### call F_6_make_function_prediction_by_MHS_and_mining_rules_and_network.pl ##
 my($F_6_MHS_net) = $dir_base."/"."F_6_make_function_prediction_by_MHS_and_mining_rules_and_network.pl";
 -s $F_6_MHS_net || die "Cannot find script $F_6_MHS_net, check the path setting\n";
 my($F_6_MHS_net_out) = $dir_output."/"."F_6_MHS_net_out";
 $return_val = system("perl $F_6_MHS_net $dir_input_BLAST 0.01 $F_4_out $mining_rule $testsequenceall $F_5_out $F_6_MHS_net_out");
 if($return_val)
 {
    print "perl $F_6_MHS_net $dir_input_BLAST 0.01 $F_4_out $mining_rule $testsequenceall $F_5_out $F_6_MHS_net_out fails!\n";
    exit(0);
 }

#### call F_6_sequence_based_prediction.pl  ###
 my($F_6_seq) = $dir_base."/"."F_6_sequence_based_prediction.pl";
 my($seq_database) = $dir_base_data."/"."new_Trained_sequence2GO";
 my($F_6_seq_out) = $dir_output."/"."F_6_seq_out";
 $return_val = system("perl $F_6_seq $seq_database $addr_fasta $F_6_seq_out");
 if($return_val)
 {
    print "perl $F_6_seq $seq_database $addr_fasta $F_6_seq_out fails!\n";
    exit(0);
 }

### call F_6_1_function_prediction_by_network_only.pl  #################
 my($F_6_net) = $dir_base."/"."F_6_1_function_prediction_by_network_only.pl";
 my($F_6_net_out) = $dir_output."/"."F_6_net_out";
 $return_val = system("perl $F_6_net $dir_input_BLAST $F_5_out $F_6_net_out");
 if($return_val)
 {
    print "perl $F_6_net $dir_input_BLAST $F_5_out $F_6_net_out fails!\n";
    exit(0);
 }   
 
#### call F_7_filter_GO_terms_for_final_prediction.pl  #################
 my($F_7_filter) = $dir_base."/"."F_7_filter_GO_terms_for_final_prediction.pl";
 my($F_6_MHS_filtered) = $dir_output."/"."Filtered_F_6_MHS_out";
 my($F_6_MHS_net_filtered) = $dir_output."/"."Filtered_F_6_MHS_net_out";
 my($F_6_seq_filtered) = $dir_output."/"."Filtered_F_6_seq_out";
 my($F_6_net_filtered) = $dir_output."/"."Filtered_F_6_net_out";
 my($go_xml) = $dir_base_data."/"."go_20130615-termdb.obo-xml";

 $return_val = system("perl $F_7_filter $F_6_MHS_out $go_xml $F_6_MHS_filtered");
 if($return_val)
 {
    print "perl $F_7_filter $F_6_MHS_out $go_xml $F_6_MHS_filtered fails!\n";
    exit(0);
 }
 $return_val = system("perl $F_7_filter $F_6_MHS_net_out $go_xml $F_6_MHS_net_filtered");
 if($return_val)
 {
    print "perl $F_7_filter $F_6_MHS_net_out $go_xml $F_6_MHS_net_filtered fails!\n";
    exit(0);
 }
 $return_val = system("perl $F_7_filter $F_6_seq_out $go_xml $F_6_seq_filtered");
 if($return_val)
 {
    print "perl $F_7_filter $F_6_seq_out $go_xml $F_6_seq_filtered fails!\n";
    exit(0);
 }

 $return_val = system("perl $F_7_filter $F_6_net_out $go_xml $F_6_net_filtered");
 if($return_val)
 {
    print "perl $F_7_filter $F_6_net_out $go_xml $F_6_net_filtered fails!\n";
    exit(0);
 }
 
#### call F_8_1_re_range_prediction.pl       ##########################
 my($F_8_rerange) = $dir_base."/"."F_8_1_re_range_prediction.pl";
 my($B_final_MHS) = $dir_output."/"."F_S3_MHS_score";
 my($B_final_MHS_net) = $dir_output."/"."F_S2_MHS_net_score";
 my($final_MHS) = $dir_output."/"."S3_MHS_score";
 my($final_MHS_net) = $dir_output."/"."S2_MHS_net_score";

 $return_val = system("perl $F_8_rerange $F_6_MHS_filtered $B_final_MHS");
 if($return_val)
 {
    print "perl $F_8_rerange $F_6_MHS_filtered $B_final_MHS fails!\n";
    exit(0);
 }
 $return_val = system("perl $F_8_rerange $F_6_MHS_net_filtered $B_final_MHS_net");
 if($return_val)
 {
    print "perl $F_8_rerange $F_6_MHS_net_filtered $B_final_MHS_net fails!\n";
    exit(0);
 }

### call F_8_2_combine_by_split.pl  ###################################
 my($F_8_2) = $dir_base."/"."F_8_2_combine_by_split.pl";
 my($F_8_2_out_dir) = $dir_output."/"."F_8_2_splitting_for_combined_score";
 $return_val = system("perl $F_8_2 $F_6_MHS_filtered 0.876 $F_6_net_filtered 0.404 $F_6_seq_filtered 0.502 $F_8_2_out_dir");
 if($return_val)
 {
    print "perl $F_8_2 $F_6_MHS_filtered 0.876 $F_6_net_filtered 0.404 $F_6_seq_filtered 0.502 $F_8_2_out_dir fails!\n";
    exit(0);
 } 

#### call F_8_3_re_range_for_combination_score.pl ####################
 my($F_8_3) = $dir_base."/"."F_8_3_re_range_for_combination_score.pl";
 my($B_final_combined) = $dir_output."/"."F_S1_MHS_net_seq_score";
 my($final_combined) = $dir_output."/"."S1_MHS_net_seq_score";
 my($F_8_2_final) = $F_8_2_out_dir."/"."Final_combined_file";
 $return_val = system("perl $F_8_3 $F_8_2_final $B_final_combined");
 if($return_val)
 {
    print "perl $F_8_3 $F_8_2_final $B_final_combined fails!\n";
    exit(0);
 }

#### add the category information for the final output ###############
 my($add_category) = $dir_base."/"."F_9_add_category.pl";
 $return_val = system("perl $add_category $B_final_combined $go_xml $final_combined");
 if($return_val)
 {
    print "perl $add_category $B_final_combined $go_xml $final_combined fails!\n";
    exit(0);
 }

 $return_val = system("perl $add_category $B_final_MHS_net $go_xml $final_MHS_net");
 if($return_val)
 {
    print "perl $add_category $B_final_MHS_net $go_xml $final_MHS_net fails!\n";
    exit(0);
 }

 $return_val = system("perl $add_category $B_final_MHS $go_xml $final_MHS");
 if($return_val)
 {
    print "perl $add_category $B_final_MHS $go_xml $final_MHS fails!\n";
    exit(0);
 }

######## draw statistical figures ###################################
 my($draw_f1) = $dir_base."/"."F_10_draw_stat_GO.pl";
 my($f1_MHS) = $dir_output."/"."S_f1_MHS.svg";
 my($f1_MHS_net) = $dir_output."/"."S_f1_MHS_net.svg";
 my($f1_MHS_net_seq) = $dir_output."/"."S_f1_MHS_net_seq.svg";
 $return_val = system("perl $draw_f1 $final_MHS $f1_MHS");
 if($return_val)
 {
    print "perl $draw_f1 $final_MHS $f1_MHS fails!\n";
    exit(0);
 }
 $return_val = system("perl $draw_f1 $final_MHS_net $f1_MHS_net");
 if($return_val)
 {
    print "perl $draw_f1 $final_MHS_net $f1_MHS_net fails!\n";
    exit(0);
 }
 $return_val = system("perl $draw_f1 $final_combined $f1_MHS_net_seq");
 if($return_val)
 {
    print "perl $draw_f1 $final_combined $f1_MHS_net_seq fails!\n";
    exit(0);
 }

 

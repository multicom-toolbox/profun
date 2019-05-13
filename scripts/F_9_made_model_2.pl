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
 sub get_specy_name($);
 if(@ARGV!=3)
 {
	 print "This script builds the model 2 for CAFA 2013 only. You can revise the key word and team name inside the scripts. MIS score is used for model2\n";
	 print "perl $0 addr_F_8_MIS_score_CAFA2013 dir_all_sequences(Downloaded from CAFA website for all targets of CAFA 2014) dir_output\n";
	 print "For example:\n";
	 print "perl $0 ../CAFA_2013_result/F_8_reranged_filtered_valid_prediction_by_MHS_mining_rules_CAFA2 ../downloaded_testing_from_CAFA/CAFA-2013-targets/all_fasta_file ../CAFA_2013_result/F_9_model2_folder\n";
         print "perl $0 ../CAFA_2013_result/F_8_reranged_filtered_valid_prediction_by_MHS_mining_rules_for_other_species_CAFA2 ../downloaded_testing_from_CAFA/CAFA-2013-other ../CAFA_2013_result/F_9_model2_for_other_species_folder\n";

	 exit(0);
 }

 my($addr_input)=$ARGV[0];
 my($dir_seq)=$ARGV[1];
 my($dir_output)=$ARGV[2];
 -s $dir_output || system("mkdir $dir_output");
 -s $dir_seq || die "cannot open fasta folder $dir_seq\n";
 my($IN,$OUT,$line,$i);
 ####### team name ##########
 my($TEAM_NAME)="profun";
 my($TEAM_ID)="profun";
 ############################
 ####### key words ##########    
 my(@key_words)=();          # comma and full stop, such as protein, gene.
 $key_words[0]="sequence-profile alignment";
 $key_words[1]="de novo prediction";


 ############################
 ###### model number ########
 my($model_number)="2";
 ############################
 
 $|=1;
 ######## process the input fasta sequences, get the target name for each species #########
 my($file,$path,$key,$value,$spe_name);

 my(%target2spe)=();        # key is target name, and value is species
 my(%all_spe)=();           # store all species names
 my(@files,@tem_split);
 opendir(DIR,"$dir_seq");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 $path=$dir_seq."/".$file;
	 print "Processing $path ...\n";
	 $spe_name=get_specy_name($file);         # get the species name from the file name
     if(not exists $all_spe{$spe_name})
	 {
		 $all_spe{$spe_name}=1;
	 }
	 $IN = new FileHandle "$path";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 if(substr($line,0,1) eq ">")
		 {# this is what we want
			 @tem_split=split(/\s+/,$line);
                         if($tem_split[0] ne ">")
                         {
			    @tem_split=split(/\>/,$tem_split[0]);
                         }
			 if(not exists $target2spe{$tem_split[1]})
			 {
				 $target2spe{$tem_split[1]}=$spe_name;
			 }
			 else
			 {
				 print "Warning, we found a duplicate $tem_split[1], check this target name! in $path\n";
			 }
		 }
	 }
	 $IN->close();
 }
 
 print "We get species name as follows:\n";
 foreach $key (keys %all_spe)
 {
	 print $key."\t";
 }
 print "\n";

 ##############################################################################################################################

 my(%check_target)=();                 # the key is the target we made predictions, and the value is the number of go terms we made for the target
 
 ##### now start processing the input file and output prediction for each prediction and also for all species ######
 foreach $key (keys %all_spe)
 {
	 print "Preparing the output model for species $key ...\n";
	 $path=$dir_output."/".$TEAM_ID."_".$model_number."_".$key.".txt";      # the output file 
	 $OUT = new FileHandle ">$path";
     ### head infor ####
	 print $OUT "AUTHOR\t$TEAM_NAME\n";
	 print $OUT "MODEL\t$model_number\n";
	 print $OUT "KEYWORDS";
	 print $OUT "\t$key_words[0]";
	 for($i=1;$i<@key_words;$i++)
	 {
		 print $OUT ", $key_words[$i]";
	 }
	 print $OUT ".\n";
     ###################
	 $IN=new FileHandle "$addr_input";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split == 4)
		 {# name target_name go score 

			 if(not exists $target2spe{$tem_split[1]})
			 {
				 print "Error, never seen this target name $tem_split[1] in the input species fasta files $addr_input, check the reason!\n";
				 exit(0);
			 }
			 if($target2spe{$tem_split[1]} ne $key)
			 {
				 next;
			 }
             ########## check target ##########
             if(not exists $check_target{$tem_split[1]})
			 {
				 $check_target{$tem_split[1]}=1;
			 }
			 else
			 {
				 $check_target{$tem_split[1]}++;
			 }
			 if($check_target{$tem_split[1]}>1500)
			 {
				 print "Warning, this target $tem_split[1], we made more than 1500 go predictions, not correct format!\n";
				 next;
			 }
             ##################################
			 print $OUT $tem_split[1]."\t".$tem_split[2]."\t".$tem_split[3]."\n"; 
		 }
		 elsif(@tem_split == 3)
		 {

			 if(not exists $target2spe{$tem_split[0]})
			 {
				 print "Error, never seen this target name $tem_split[0] in the input species fasta files $addr_input, check the reason!\n";
				 exit(0);
			 }
			 if($target2spe{$tem_split[0]} ne $key)
			 {
				 next;
			 }
             ########## check target ##########
             if(not exists $check_target{$tem_split[0]})
			 {
				 $check_target{$tem_split[0]}=1;
			 }
			 else
			 {
				 $check_target{$tem_split[0]}++;
			 }
			 if($check_target{$tem_split[0]}>1500)
			 {
				 print "Warning, this target $tem_split[0], we made more than 1500 go predictions, not correct format!\n";
				 next;
			 }
             ##################################
			 print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\n"; 
		 }
		 else
		 {
			 print "Warning !!!!!!! skip $line\n";
		 }
	 }
	 $IN->close();

	 ### end infor ###
	 print $OUT "END";
	 #################
	 $OUT->close();
 }



 ############ now in case we needed, we also made one model with all species prediction #########



 $path=$dir_output."/".$TEAM_ID."_all_species_prediction".".txt";      # the output file 
 $OUT = new FileHandle ">$path";
 ### head infor ####
 print $OUT "AUTHOR\t$TEAM_NAME\n";
 print $OUT "MODEL\t$model_number\n";
 print $OUT "KEYWORDS";
 print $OUT "\t$key_words[0]";
 for($i=1;$i<@key_words;$i++)
 {
	 print $OUT ", $key_words[$i]";
 }
 print $OUT ".\n";
 ###################
 $IN=new FileHandle "$addr_input";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 @tem_split=split(/\s+/,$line);
	 if(@tem_split == 4)
	 {# name target_name go score 
		 if(not exists $target2spe{$tem_split[1]})
		 {
			 print "Error, never seen this target name $tem_split[1] in the input species fasta files $addr_input, check the reason!\n";
			 exit(0);
		 }
		 print $OUT $tem_split[1]."\t".$tem_split[2]."\t".$tem_split[3]."\n"; 
	 }
	 elsif(@tem_split == 3)
	 {
		 if(not exists $target2spe{$tem_split[0]})
		 {
			 print "Error, never seen this target name $tem_split[0] in the input species fasta files $addr_input, check the reason!\n";
			 exit(0);
		 }
		 print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\n"; 
	 }
	 else
	 {
		 print "Warning !!!!!!! skip $line\n";
	 }
 }
 $IN->close();

 ### end infor ###
 print $OUT "END";
 #################
 $OUT->close();

 ########## check the missing target ###########
 foreach $key (keys %target2spe)
 {
	 if(not exists $check_target{$key})
	 {
		 print "Missing prediction for target $key!\n";
	 }
 }
 


 sub get_specy_name($)
 {
	 my($file_name)=@_;
	 my(@tem)=split(/\./,$file_name);
	 if(@tem != 3)
	 {
		 print "Check the input file name $file_name, is not the correct format for CAFA 2013 release fasta sequence file name!\n";
		 $tem[1]="NULL";
	 }
	 my($s)=$tem[1];
	 return $s;
 }

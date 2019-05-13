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
 use Cwd;
 use Cwd 'abs_path';
 sub get_frequency($);

 use POSIX qw/ceil/; #this is for function ceil()
 if (@ARGV != 10)
 { # @ARGV used in scalar context = number of args

      print "This program do analysis on the aligned true and predicted gene function, use the tool fastsemsim to evaluate the similarity, filtered the result and convert similarity score into bins, use R draw pictures.!\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 6/30/2013\n\n");
	  print("***********************************************************\n");
	  print "parameters:\n1.aligned_true_and_prediction(in this folder, each file is one alignment) \n";
	  print "2. Dir of fastsemsim folder. (We know the script name is ./fastSemSim.sh)\n";
	  print "3. Address of script to filter fastsemsim result. \n";
	  print "4. Address of script to calculate the similarity between two gene's go function.\n";
	  print "5. Address of script to get equal width of similarity.\n";
	  print "6. Dir of fastsemsim processed folder\n";
	  print "7. Dir of filtered fastsemsim result\n";
	  print "8. Dir of divide bin\n";
	  print "9. Dir of binned score of similarity score\n";
	  print "10. Dir of binned score of similarity in pdf \n";

	  print "\n\nFor example:\n";

          print "perl $0 ../result/5_aligned_true_and_predictions /exports/store1/tool/fastSemSim-0.7 /exports/store1/rcrg4/special_project/4_10_script/0_2_filter_fastsemsim_result.pl /exports/store1/rcrg4/special_project/4_10_script/6_get_similarity_of_GO_of_gene_pairs.pl /exports/store1/rcrg4/special_project/4_10_script/7_get_equal_width_distance_for_go.pl ../result/6_one_step_analysis/fastsemsim_processed_result ../result/6_one_step_analysis/filtered_fastsemsim_result ../result/6_one_step_analysis/divide_bin_result ../result/6_one_step_analysis/bin_result ../result/6_one_step_analysis/bin_pdf_result\n";


	  exit(0);
 }
 my($input_align)=abs_path($ARGV[0]);
 my($dir_tool)=$ARGV[1];
 my($script_filter)=abs_path($ARGV[2]);
 my($script_cal)=abs_path($ARGV[3]);
 my($script_equal)=abs_path($ARGV[4]);
 my($dir_tool_pro)=abs_path($ARGV[5]);
 my($dir_tool_filtered)=abs_path($ARGV[6]);
 my($dir_divide_bin)=abs_path($ARGV[7]);
 my($dir_bin)=abs_path($ARGV[8]);
 my($dir_pdf)=abs_path($ARGV[9]);

 -s $input_align || die "cannot open input true and prediction aligned file:  $input_align\n";
 -s $dir_tool || die "cannot open input fastsemsim folder $dir_tool\n";
 -s $script_filter || die "cannot open input $script_filter\n";
 -s $script_cal || die "cannot open input $script_cal\n";
 -s $dir_tool_pro || system("mkdir $dir_tool_pro");
 -s $dir_tool_filtered || system("mkdir $dir_tool_filtered");
 -s $dir_divide_bin || system("mkdir $dir_divide_bin");
 -s $dir_bin || system("mkdir $dir_bin");
 -s $dir_pdf || system("mkdir $dir_pdf");

 my($target,$key,$value,$path_input,$path_output,$IN,$OUT,$line,$i,$j,$path_real_tool,$return_val);
 my(@tem_split,@targets,@tem_2);
 
 $path_real_tool=$dir_tool."/fastSemSim.sh";           # the real path to run fastsemsim 
 
 print "Use fastsemsim to process data ...\n";
 ####### 1. Use fastsemsim to processed the aligned data #########
 opendir(DIR,"$input_align");
 @targets=readdir(DIR);
 foreach $target (@targets)
 {
	 if($target eq '.' || $target eq '..')
	 {
		 next;
	 }
	 $path_input=$input_align."/".$target;
	 chdir($dir_tool);       # go to the tool folder
     $path_output=$dir_tool_pro."/".$target;
	 $return_val=system("$path_real_tool --GOTerms -s G-SESAME -q $path_input --pairs --consider_has_part --go /exports/store1/tool/fastSemSim-0.7/examples/go_daily-termdb.obo-xml.gz > $path_output");
	 if($return_val !=0)
	 {
		 print "$path_real_tool --GOTerms -s G-SESAME -q $path_input --pairs --consider_has_part --go /exports/store1/tool/fastSemSim-0.7/examples/go_daily-termdb.obo-xml.gz > $path_output fails!\n";
		 next;
	 }
 }
 print "Filter data ...\n";
 ##### 2. Filtered fastsemsim processed result, remove the go pairs with no similarity, not in the same category  ######
 opendir(DIR,"$dir_tool_pro");
 @targets=readdir(DIR);
 foreach $target (@targets)
 {
	 if($target eq '.' || $target eq '..')
	 {
		 next;
	 }
	 $path_input=$dir_tool_pro."/".$target;
	 $path_output=$dir_tool_filtered."/".$target;
	 $return_val=system("perl $script_filter $path_input $path_output");
	 if($return_val !=0)
	 {
		 print "perl $script_filter $path_input $path_output fails!\n";
		 next;
	 }
 }
 
 print "Calculate similarity ...\n";
 ##### 3. calculate the gene gene similarity ########
 opendir(DIR,"$dir_tool_filtered");
 @targets=readdir(DIR);
 foreach $target (@targets)
 {
	 if($target eq '.' || $target eq '..')
	 {
		 next;
	 }
	 ##### first use ave metric #######
	 $path_input=$dir_tool_filtered."/".$target;
	 $path_output=$dir_divide_bin."/".$target."_ave";
	 $return_val=system("perl $script_cal $path_input ave $path_output");
	 if($return_val !=0)
	 {
		 print "perl $script_cal $path_input ave $path_output fails!\n";
		 next;
	 }
	 ##### use max metric #######
	 $path_input=$dir_tool_filtered."/".$target;
	 $path_output=$dir_divide_bin."/".$target."_max";
	 $return_val=system("perl $script_cal $path_input max $path_output");
	 if($return_val !=0)
	 {
		 print "perl $script_cal $path_input max $path_output fails!\n";
		 next;
	 }
 }
 print "Get equal width bin score ...\n";
 ##### 4. Get equal width bin score ########
 opendir(DIR,"$dir_divide_bin");
 @targets=readdir(DIR);
 foreach $target (@targets)
 {
	 if($target eq '.' || $target eq '..')
	 {
		 next;
	 }
	 $path_input=$dir_divide_bin."/".$target;
	 $path_output=$dir_bin."/".$target;
	 $return_val=system("perl $script_equal $path_input 0.1 $path_output");
	 if($return_val !=0)
	 {
		 print "perl $script_equal $path_input 0.1 $path_output fails!\n";
		 next;
	 }
 }


 print "Use R draw pictures ...\n";


 ##### 5. Read the binned score, and use R to draw the pdf file ######
 opendir(DIR,"$dir_bin");
 my($path_tmp);
 my(@frequencies);
 @targets=readdir(DIR);
 foreach $target (@targets)
 {
	 if($target eq '.' || $target eq '..')
	 {
		 next;
	 }
	 $path_input=$dir_bin."/".$target;      
     $path_output=$dir_pdf."/".$target.".pdf";

	 @frequencies = ();
	 @frequencies=get_frequency($path_input);       # from the binned similarity score input file, get the frequency for each bin
     ########## write the R script to draw pictures #########
     $path_tmp=$dir_pdf."/".$target.".r";
	 $OUT=new FileHandle ">$path_tmp";
	 print $OUT "y<-c($frequencies[0]";
	 for($i=1;$i<@frequencies;$i++)
	 {
		 print $OUT ",".$frequencies[$i];
	 }
	 print $OUT ")\n";
	 print $OUT "pdf(\"$path_output\")\n";
	 print $OUT "barplot(y,xlab=\"similarity score in bin\",ylab=\"Gene pair frequency\")\n";
	 print $OUT "dev.off()\n";
	 $OUT->close();
	 system("R --vanilla < $path_tmp");
 }
 
sub get_frequency($)
{
	my($input)=@_;
	my($IN,$line);
	my(@scores)=();
	my($index)=0;
	my(@tem);
	$IN = new FileHandle "$input";
	defined($IN) || die "Cannot open $input\n";
    if(defined($line=<$IN>))
	{# skip the first line
	}
	while(defined($line=<$IN>))
	{
		chomp($line);
		$line=~s/\s+$//;
		@tem=split(/\s+/,$line);
		if(@tem<2)
		{
			next;
		}
		$scores[$index++]=$tem[1];
	}
	$IN->close();
	return @scores;
}

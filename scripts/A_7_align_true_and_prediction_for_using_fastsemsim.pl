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

      print "This program aligns the true and prediction of gene functions, prepare for fastsemsim to make prediction!\n";
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 6/27/2013\n\n");
	  print("***********************************************************\n");
	  print "perl $0 real_gene_function predicted_gene_function addr_output\n";
	  print "For example:\n";
	  print "perl $0 ../result/2_real_go_function_for_each_protein/other.mitab.go_functions ../result/4_protein_function_predictions_from_one_step_neighbours/other.function_prediction ../result/5_aligned_true_and_predictions/aligned_for_fastsemsim.others\n";
          exit(0);
 }
 my($input_real)=$ARGV[0];
 my($input_predict)=$ARGV[1];
 my($output)=$ARGV[2];
 -s $input_real || die "cannot open $input_real\n";
 -s $input_real || die "cannot open $input_predict\n";

 my($IN,$OUT,$line);
 my(@tem_split);
 
 ######### 1 read the input real gene functions #########
 my(%hash_real)=();       
 $IN=new FileHandle "$input_real";
 defined($IN) || die "Cannot open input $input_real!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\|/,$line);
	 if(@tem_split<2)
	 {
		 next;
	 }
	 if(not exists $hash_real{$tem_split[0]})
	 {
		 $hash_real{$tem_split[0]}=$tem_split[1];
	 }
	 else
	 {
		 print "Why $tem_split[0] duplicate??/ check $line\n";
	 }
 }
 $IN->close();
 
 ######### 2. read predicted gene function and print out ##########
 $OUT=new FileHandle ">$output";
 
 my(@real_go,@pre_go);
 my($i,$j);
 $IN=new FileHandle "$input_predict";
 defined($IN) || die "Cannot open input $input_predict!\n";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\|/,$line);
	 if(@tem_split<2)
	 {
		 next;
	 }
	 if(not exists $hash_real{$tem_split[0]})
	 {
		 print "Skip gene $tem_split[0], since we don't get the real gene function!\n";
		 next;
	 }
	 ####### get the real go and predict go ######
	 @real_go=split(/\s+/,$hash_real{$tem_split[0]});
	 @pre_go=split(/\s+/,$tem_split[1]);

#print "@real_go\n";
#print "@pre_go\n";


     for($i=0;$i<@real_go;$i++)
	 {
		 for($j=0;$j<@pre_go;$j++)
		 {
			 print $OUT $real_go[$i]."\t".$pre_go[$j]."\n";
#print $real_go[$i]."\t".$pre_go[$j]."\n";

		 }
	 }
	 print $OUT "#ENDOFGENE\t#ENDOFGENE\ttag\n";

#exit(0);

 }
 $IN->close();
 $OUT->close();
 

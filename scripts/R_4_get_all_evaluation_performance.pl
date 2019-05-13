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
 use POSIX qw/ceil/;

 if (@ARGV != 2)
 {
         print "This program collects all the evaluation result by R_3 ..., and get the final average precision and recall for 5 cross validation. \n";
         print "We should also consider the number of rules, the computation complexity! Finally, s=0.05, c=90 is ok for us!\n";
         print "0.02 => 171817 rules, 0.03 => 120114 rules, 0.04 => 62707, 0.05 => 51356 trainning rules!\n";
 
		 print "perl $PROGRAM_NAME dir_input(each folder has file Performance inside store the precision and recall) addr_output\n";
		 print "For example:\n";
		 print "perl $PROGRAM_NAME ../result/R_3_evaluation_results ../result/R_4_overall_performance\n";


	 exit(0);
 }
 my($dir_input)=$ARGV[0];
 my($output)=$ARGV[1];

 my($IN,$OUT,$line,$file,$path_folder,$name,$key,$value,$i);
 my(@files,@tem_split);

 my(%hash_cross)=();           # key is the target_name of cross validation, value is the path for each cross validation

 opendir(DIR,"$dir_input");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 $path_folder=$dir_input."/".$file."/"."Performance";
	 @tem_split=split(/\_/,$file);
	 if(@tem_split<2)
	 {
		 print "The name is not correct, maybe ,check $dir_input/$file\n";
		 exit(0);
	 }
	 $name=$tem_split[0];
     for($i=1;$i<@tem_split-2;$i++)
	 {
		 $name.="_".$tem_split[$i];
	 }
	 if(exists $hash_cross{$name})
	 {
		 $hash_cross{$name}.="|".$path_folder;
	 }
	 else
	 {
		 $hash_cross{$name}=$path_folder;
	 }
 }

 $OUT = new FileHandle ">$output";
 my(@tem2);
 my(%hash_score)=();

 my($ave_recall,$index,$ave_precision,$tmp_rec,$tmp_pre);
 foreach $key (keys %hash_cross) 
 {
	 @tem_split=split(/\|/,$hash_cross{$key});
	 $ave_recall=0;
	 $ave_precision=0;
	 $index=0;

	 for($i=0;$i<@tem_split;$i++)
	 {
		 $tmp_pre=-1;
		 $tmp_rec=-1;
		 $index++;
		 #print $tem_split[$i]."\n";
		 if(!-s $tem_split[$i])
		 {
			 print "Not existing $tem_split[$i]!\n";
			 next;
		 }
		 $IN=new FileHandle "$tem_split[$i]";
		 while(defined($line=<$IN>))
		 {
			 chomp($line);
			 if($line =~ m/precision/)
			 {
				 @tem2=split(/\s+/,$line);
                 $tmp_pre=$tem2[2];
				 $ave_precision+=$tmp_pre;
			 }
			 elsif($line =~ m/recall/)
			 {
				 
				 @tem2=split(/\s+/,$line);
				 $tmp_rec=$tem2[2];
				 $ave_recall+=$tmp_rec;
			 }
		 }
		 $IN->close();
	     if($tmp_pre == -1 || $tmp_rec == -1)
		 {
			 print "Error, the precision or recall is not exists in $tem_split[$i]!\n";
			 exit(0);
		 }
		 #print $OUT "$tem_split[$i](precision = $tmp_pre, recall = $tmp_rec);";
	 }
	 $ave_precision=$ave_precision/$index;
	 $ave_recall=$ave_recall/$index;
	 print $OUT "$key ";
	 print $OUT " ==> average precision = $ave_precision, average recall = $ave_recall. Multiplication = ";
	 $tmp_pre=$ave_precision*$ave_recall;
	 print $OUT $tmp_pre."\n";

	 if(not exists $hash_score{$key})
	 {
		 $hash_score{$key} = $ave_precision."|".$ave_recall."|".$tmp_pre;
	 }

 }
 print $OUT "\n\n#### support, confidence, precision, recall, multiplication of precision and recall ################################\n";
 foreach $key (keys %hash_score)
 {
	 @tem_split=split(/\_/,$key);
	 @tem2=split(/\|/,$hash_score{$key});
	 print $OUT $tem_split[1]."\t".$tem_split[3]."\t".$tem2[0]."\t".$tem2[1]."\t".$tem2[2]."\n";
 }

 $OUT->close();


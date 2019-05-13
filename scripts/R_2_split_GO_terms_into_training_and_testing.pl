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
 use List::Util qw/shuffle/;
 sub get_num_line($);

 if (@ARGV != 3)
 {
         print "This program will create a ten(or other) cross validation from a input dataset\n";
		 print "perl $PROGRAM_NAME addr_input(uniprot_parsed_GO_terms) num_of_cross_validation dir_output \n";
		 print "For example:\n";
		 print "perl $PROGRAM_NAME ../result/R_1_mining_rules_from_uniprot/uniprot_parsed_GO_terms 5 ../result/R_2_cross_5_validation_for_GO_terms\n";

	 exit(0);
 }
 my($addr_input)=$ARGV[0];
 my($num)=$ARGV[1];
 my($dir_output_SVM)=$ARGV[2];

 -s $dir_output_SVM || system("mkdir $dir_output_SVM");

 ###################################################################
 my($line,$IN,$OUT,$path,$file,$OUT2,$OUT3,$format_line);
 my(@tem_split,@files);
 my($count_traing_data);
 -s $addr_input || die "Cannot open $addr_input\n";
 -s $dir_output_SVM || system("mkdir $dir_output_SVM");

 if($num<0)
 {
	 print "The input number $num should be larger than 0, this is the number of cross validation you need, usually we use 10 or 5\n";
	 exit(0);
 }
     my($total)=0;

	 ######## read the input file ###########
	 $IN=new FileHandle "$addr_input";
	 defined($IN) || die "Cannot open input $addr_input\n";
     while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split<1)
		 {
			 next;
		 }
		 $total++;

	 }
     $IN->close();
	 my($interval)=ceil($total/$num);
  print "There are $total number of lines in the input file $addr_input, we do $num cross-validation now, each part has $interval lines ...\n";
     
     my(@ranks);
	 my(@shuffle_target);
	 my($i,$j);
     @ranks=();
     for($i=0;$i<$total;$i++)
     {
	    $ranks[$i]=$i;
     }
     @shuffle_target=();
     @shuffle_target = (shuffle(@ranks))[0..$total-1];


	 my($start)=0;
	 my($cross_num)=0;     
	 my($line_num,$end2);
     my(%hash_detect)=(); 
  
     my($path_train,$path_test,$path_train2,$path_test2,$path_train3,$path_test3);
	 for($cross_num=0;$cross_num<$num;$cross_num++)
	 {######## start from $cross_num * interval, to ($cross_num+1) * interval as testing data, others are training data 
        $path_train=$dir_output_SVM."/"."training_".$cross_num;   # training dataset SVM format
		$path_test=$dir_output_SVM."/"."testing_".$cross_num;     # testing dataset  SVM format
        


		$start=$cross_num*$interval;                          # this is the starting point
		while($start>=$total)
		{
			$start=$start-$total;
		}
		%hash_detect=();



		for($j=0;$j<$interval;$j++)
		{# get $interval number of lines for the testing dataset
			$line_num=$start+$j;                # this is the line number
		    while($line_num>=$total)
		    {
			    $line_num=$line_num-$total;
		    }
			if(not exists $hash_detect{$shuffle_target[$line_num]})
			{
				$hash_detect{$shuffle_target[$line_num]}=1;
			}
		}
		######### output the testing dataset ###########
		$OUT=new FileHandle ">$path_test";        # for SVM


		$IN=new FileHandle "$addr_input";
		$i=0;
	    defined($IN) || die "Cannot open input $addr_input\n";
        while(defined($line=<$IN>))
	    {
		   chomp($line);
		   $line=~s/\s+$//;
		   @tem_split=split(/\s+/,$line);
		   if(@tem_split<1)
		   {
			 next;
		   }
		   if(exists $hash_detect{$i})
		   {
			   print $OUT $line."\n";
		   }
           $i++;
		}
		$IN->close();
		$OUT->close();
		$count_traing_data=0;
		#################################################
		########## output the training dataset ##########
		%hash_detect=();
		print "Now the line_num is $line_num, the start position is $start, the total is $total\n";
		for(;$line_num!=$start;$line_num++)
		{
		    while($line_num>=$total)
		    {
			    $line_num=$line_num-$total;
				$line_num--;                              # since the for(;$line_num!=$start;$line_num++) will increase $line_num automaticly, so here we need to decrease 1
		    }
			if(not exists $hash_detect{$shuffle_target[$line_num]})
			{
				$hash_detect{$shuffle_target[$line_num]}=1;
			}
		    $count_traing_data++;
		}

	
	    ######### output the testing dataset ###########
		$OUT=new FileHandle ">$path_train";        # for SVM


		$IN=new FileHandle "$addr_input";
		$i=0;
	    defined($IN) || die "Cannot open input $addr_input\n";
        while(defined($line=<$IN>))
	    {
		   chomp($line);
		   $line=~s/\s+$//;
		   @tem_split=split(/\s+/,$line);
		   if(@tem_split<1)
		   {
			 next;
		   }
		   if(exists $hash_detect{$i})
		   {
			   print $OUT $line."\n";

		   }
           $i++;
		}
		$IN->close();
		$OUT->close();
		#################################################

	 }


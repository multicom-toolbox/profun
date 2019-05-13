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



  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
          print "Not used, since the performance is not very good!\n";
	  print "This script will re-range the prediction data, score will start from the maximum score of that target, decrease by 0.01 each time. We want to get better threshold metric evaluation result! \n";
	  print "perl $0 addr_input addr_output\n";
	  print "For example:\n";
	  print "perl $0 ../result/processed_prediction_F_6_MHS_mining_rule_CAFA2 ../result/From_max_to_re_ranged_prediction_by_MHS_and_mining_rule_CAFA2\n";

	  exit(0);
	}

  my($input)=$ARGV[0];
  my($output)=$ARGV[1];

  my($IN,$OUT,$line,$key,$value,$i,$pre_6_deci,$new_6_deci);
  my(@tem_split,@tem2);
  my(%hash)=();                          # key is target name, value is modelname targetname GO currentscore|mappedscore
  $OUT = new FileHandle ">$output";
  $IN=new FileHandle "$input";
  while(defined($line=<$IN>))
  {
	  chomp($line);
	  @tem_split=split(/\s+/,$line);

	  if(@tem_split == 4)
	  {#modelname targetname GO score
		  $new_6_deci=$tem_split[3];                              # new get 6 decimal value 

          $tem_split[3]=sprintf("%.2f",$tem_split[3]);            # convert to 0.01 - 0.99 format
		  if(not exists $hash{$tem_split[1]})
		  {
			  $hash{$tem_split[1]} = $line."|".$tem_split[3];
		  }
		  else
		  {# already have some function
			  @tem2=split(/\|/,$hash{$tem_split[1]});
			  @tem2=split(/\s+/,$tem2[0]);
			  $pre_6_deci=$tem2[3];          # we store the 6 decimal value previously generated

			  $tem2[3]=sprintf("%.2f",$tem2[3]);
			  if($pre_6_deci < $new_6_deci)
			  {
				  print "Error, the data is not ranked, check $line, $pre_6_deci and $new_6_deci!\n";
				  exit(0);
			  }
#			  if($tem_split[3]>$tem2[3])
#			  {
#				  print "The input prediction is not ranked! Check here! $line, previous get $tem2[3], now is $tem_split[3]\n";
#				  exit(0);
#			  }
			  elsif($new_6_deci<$pre_6_deci)
			  {# the new is less
				  @tem2=split(/\|/,$hash{$tem_split[1]});
				  $tem2[1]-=0.01;
				  if($tem2[1]<=0.01)
				  {
					  print "Warning,  check why it decrease to negative??? $line,  We set it to minimum value 0.01\n";
					  $tem2[1]=0.01;
					  #exit(0);
				  }
				  $hash{$tem_split[1]}=$line."|".$tem2[1];
			  }
			  else
			  {# equal
				  @tem2=split(/\|/,$hash{$tem_split[1]});
				  $hash{$tem_split[1]}=$line."|".$tem2[1];
			  }
		  }
		  @tem2=split(/\s+/,$line);
		  print $OUT $tem2[0]."\t".$tem2[1]."\t".$tem2[2]."\t";
		  @tem2=split(/\|/,$hash{$tem_split[1]});
          $tem2[1]=sprintf("%.2f",$tem2[1]);
		  print $OUT $tem2[1]."\n";

	  }
	  elsif(@tem_split == 3)
	  {# targetname GO score
		  $new_6_deci=$tem_split[2];                              # new get 6 decimal value 
		  $tem_split[2]=sprintf("%.2f",$tem_split[2]);            # convert to 0.01 - 0.99 format
		  if(not exists $hash{$tem_split[0]})
		  {
			  $hash{$tem_split[0]} = $line."|".$tem_split[2];
		  }
		  else
		  {# already have some function
			  @tem2=split(/\|/,$hash{$tem_split[0]});
			  @tem2=split(/\s+/,$tem2[0]);
			  $pre_6_deci=$tem2[2];          # we store the 6 decimal value previously generated
			  $tem2[2]=sprintf("%.2f",$tem2[2]);
			  if($pre_6_deci < $new_6_deci)
			  {
				  print "The input prediction is not ranked! Check here! $line\n";
				  exit(0);
			  }
			  elsif($new_6_deci<$pre_6_deci)
			  {
				  @tem2=split(/\|/,$hash{$tem_split[0]});
				  $tem2[1]-=0.01;
				  if($tem2[1]<=0.01)
				  {
					  print "warning, check why it decrease to negative??? $line, $pre_6_deci, $new_6_deci, we set it to 0.01\n";
					  $tem2[1]=0.01;
				  }
				  $hash{$tem_split[0]}=$line."|".$tem2[1];
			  }
			  else
			  {
				  @tem2=split(/\|/,$hash{$tem_split[0]});
				  $hash{$tem_split[0]}=$line."|".$tem2[1];
			  }
		  }
		  @tem2=split(/\s+/,$line);
		  print $OUT $tem2[0]."\t".$tem2[1]."\t";
		  @tem2=split(/\|/,$hash{$tem_split[0]});
		  $tem2[1]=sprintf("%.2f",$tem2[1]);
		  print $OUT $tem2[1]."\n";
	  }
  }
  $IN->close();
  $OUT->close();

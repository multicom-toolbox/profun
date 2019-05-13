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
      print "This script is the final validation for the predictions, the input can be a model !\n";
	  print "This program will filter the GO terms for the final function prediction. We will load the valid GO terms from the file go_20130615-termdb.obo-xml, which is the version of 06/15/2013. For all the GO terms which is not valid, we will delete them!\n";
	  print "For each target, we only keep maximum 1500 GO predictions, and filter all others! The input GO predictions should be ranked from high score to low before using this script\n";
	  print "perl $0 addr_input addr_GO_term_database addr_output\n";
	  print "for example:\n";
      print "perl $0 ../test/donet_dir/final_model_2_2014.txt ../data/go_20130615-termdb.obo-xml ../test/donet_dir/processed_final_model2_2014.txt\n";
	  
	  exit(0);
	}
   my($input)=$ARGV[0];
   my($addr_GO)=$ARGV[1];
   my($output)=$ARGV[2];

   $|=1;
   my($IN,$OUT,$line,$name,$term_start,$is_obsolete,$name_space,$GO_id,$key);
   my(@tem_split);
   my(%valid_GO)=();           # get all valid GO terms .
   
   $term_start=0;
   $is_obsolete=0;
   $name_space="NULL";
   $GO_id="NULL";
   print "loading all valid GO terms ...";
##############################################################################
#
#
#               1. Loading all valid GO terms
#
#
##############################################################################
   $IN=new FileHandle "$addr_GO";
   while(defined($line=<$IN>))
   {
	   chomp($line);
	   @tem_split=split(/\</,$line);
	   @tem_split=split(/\>/,$tem_split[1]);
	   $name=$tem_split[0];           # this is the name ID 
	   if($name eq "term")
	   {# the new term start
		   $term_start=1;
	   }
	   elsif($name eq "id")
	   {# get the id
		   @tem_split=split(/\>/,$line);
		   @tem_split=split(/\</,$tem_split[1]);
		   $GO_id=$tem_split[0];
	   }
	   elsif($name eq "is_obsolete")
	   {
		   $is_obsolete=1;
	   }
	   elsif($name eq "namespace")
	   {# get the namespace
		   @tem_split=split(/\>/,$line);
		   @tem_split=split(/\</,$tem_split[1]);
		   $name_space=$tem_split[0];
	   }
	   elsif($name eq "/term")
	   {# this is end of the term
           if($name_space eq "NULL" || $GO_id eq "NULL")
		   {
			   print "Surprised, check $GO_id, $name_space, empty???\n";
			   exit(0);
		   }
		   if($term_start!=1)
		   {
			   print "Surprised, this term $GO_id, $name_space, without starting, come to the end!!!\n";
			   exit(0);
		   }
		   if($is_obsolete == 0)
		   {# now we add this valid GO term
			   if(not exists $valid_GO{$GO_id})
			   {
				   $valid_GO{$GO_id}=$name_space;
			   }
			   else
			   {
				   print "This Go term $GO_id exists more than 1 times!\n";
			   }
		   }
		   $is_obsolete=0;
           $term_start=0;
           $name_space="NULL";
           $GO_id="NULL";
	   }

   }
   $IN->close();

#   foreach $key (keys %valid_GO)
#   {
#	   print $key."\t".$valid_GO{$key}."\n";
#   }
   print "loading all predictions and filter the predictions ...";
##############################################################################
#
#
#               2. Loading input predictions
#
#
##############################################################################
  my(%target)=();             # each target should have maximum 1500 predictions
  my(@tem_2);
  $OUT = new FileHandle ">$output";
  defined($OUT) || die "Cannot open $output\n";
   $IN=new FileHandle "$input";
   while(defined($line=<$IN>))
   {
	   chomp($line);
	   @tem_split=split(/\s+/,$line);
       ######## skip the head infor ########
	   if($line =~ m/AUTHOR/ || $line =~ m/MODEL/ ||$line =~ m/KEYWORDS/ || $line =~ m/END/)
	   {
		   print "Get $line, pass by !\n";
		   print $OUT $line."\n";
		   next;
	   }

	   if(@tem_split<2)
	   {
		   print "Warning !!! skip $line \n";
		   next;
	   }
	   if(@tem_split>3)
	   {# we have four columns case, Sequence_based  T992870000001   GO:0019521      1.00000
		   $tem_split[0]=$tem_split[1];
		   $tem_split[1]=$tem_split[2];
		   $tem_split[2]=$tem_split[3];
	   }
       if($tem_split[2] == 0)
	   {
		   print "$line, Severe warning, the score should not be 0!\n";
		   exit(0);
	   }


	   if(not exists $valid_GO{$tem_split[1]})
	   {# skip the not valid GO terms
		   print "For GO terms, Skip $tem_split[1] \n";
		   next;
	   }
       if(not exists $target{$tem_split[0]})
	   {
		   $target{$tem_split[0]}="1|$tem_split[2]";
	   }
	   else
	   {
		   @tem_2=split(/\|/,$target{$tem_split[0]});
		   if($tem_2[0]>=1500)
		   {
			   print "Warning !!!!!!  Exceed the maximum number 1500 GO terms, skip rest! $line\n";
			   next;
		   }
		   $tem_2[0]++;
		   if($tem_split[2]>$tem_2[1])
		   {
			   print "Wrong, the input should ranked by the score!, check $line!\n";
			   exit(0);
		   }
		   $target{$tem_split[0]}=$tem_2[0]."|".$tem_split[2];
		   
	   }
	   print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\n";

       
   }
   $IN->close();
  $OUT->close();

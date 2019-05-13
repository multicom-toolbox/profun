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
      print "This program will filter the GO terms for the final function prediction. We will load the valid GO terms from the file go_20130615-termdb.obo-xml, which is the version of 06/15/2013. We will add the category information at the last column!\n";
	  print "perl $0 addr_input addr_GO_term_database addr_output\n";
	  print "for example:\n";
      
          print "perl $0 ../CAFA_2013_result/F_6_prediction_by_MHS_and_mining_and_network_rules_for_other_species_CAFA2 ../data/go_20130615-termdb.obo-xml ../CAFA_2013_result/F_7_filtered_valid_prediction_by_MHS_and_mining_rule_and_network_for_other_species_CAFA2\n";
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
	   if(@tem_split<2)
	   {
		   print "skip $line \n";
		   next;
	   }
	   if(@tem_split>3)
	   {# we have four columns case, Sequence_based  T992870000001   GO:0019521      1.00000
		   $tem_split[0]=$tem_split[1];
		   $tem_split[1]=$tem_split[2];
		   $tem_split[2]=$tem_split[3];
	   }
	   if(not exists $valid_GO{$tem_split[1]})
	   {# skip the not valid GO terms
		   #print "Skip $tem_split[1] \n";
		   next;
	   }
	   print $OUT $tem_split[0]."\t".$tem_split[1]."\t".$tem_split[2]."\t".$valid_GO{$tem_split[1]}."\n";

       
   }
   $IN->close();
  $OUT->close();

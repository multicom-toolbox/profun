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
      print("This program will extract the protein AC ID from the swiss-prot BLAST result, the file to use is *.list. Format is sp|P19668|BGAL_BACST Beta-galactosidase 1 OS=Bacillus stearother...   976   0.0, ...\n");
      print("***********************************************************\n\n");
	  print("Renzhi Cao, 10/5/2013\n\n");
	  print("***********************************************************\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_BLAST_result e_value_threshold addr_output!\n");
      print("\n\n***********************************************************\n");
	  print("Examples:\n");
      print("perl $PROGRAM_NAME ../CAFA_2012_prediction 0.01 ../result/F_1_all_proteins_from_BLAST\n\n");

	  exit(1) ;
    }
  my $starttime = localtime();
  print "\nThe time started at : $starttime.\n";

 my($dir_input)=$ARGV[0];
 my($threshold)=$ARGV[1];
 my($addr_output)=$ARGV[2];

 -s $dir_input || die "No input file $dir_input!\n";
 my($key,$value,$IN,$OUT,$path_input,$target_name,$path_output,$line,$file,$i,$j,$evalue);
 my(@tem_split,@tem_2,@files);
 


 my(%hash_protein)=();                      # store all protein ID , key is protein ID, value is evalue1_targetname1|evalue2_target_name2 ...
 opendir(DIR,"$dir_input");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 @tem_split=split(/\./,$file);
	 $target_name=$tem_split[0];
	 if(@tem_split<2 || $tem_split[1] ne "list")
	 {
		 print "the file name is not *.list format, skip $file !\n";
		 next;
	 }
	 $path_input=$dir_input."/".$file;

print "PRocssing $path_input !\n";
	 
	 $IN=new FileHandle "$path_input";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;

		 @tem_split=split(/\s+/,$line);
		 if(@tem_split<3)
		 {
			 next;
		 }
         $evalue=$tem_split[@tem_split-1];
		 if(substr($evalue,0,1) eq "e" || substr($evalue,0,1) eq "E")
		 {
			 $evalue="1".$evalue;
		 }

         if($evalue>$threshold)
		 {#  the evalue is too large, we don't want this hit
			 next;
		 }

		 @tem_split=split(/\|/,$line);
		 $key=$tem_split[1];
		 if(not exists $hash_protein{$key})
		 {# this is the first time hit protein
			 $hash_protein{$key}=$evalue."_".$target_name;
		 }
		 else
		 {# this is more than first time hit, in different species
			 $hash_protein{$key}.="|".$evalue."_".$target_name;
		 }
		 
	 }
	 $IN->close();
 }
 $OUT=new FileHandle ">$addr_output";
 foreach $key (keys %hash_protein)
 {
	 print $OUT $key."\t".$hash_protein{$key}."\n";
 }
 $OUT->close();
 

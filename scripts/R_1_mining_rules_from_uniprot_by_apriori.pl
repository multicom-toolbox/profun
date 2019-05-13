 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 sub filter_file($$);
 use POSIX qw/ceil/; #this is for function ceil()
   if (@ARGV != 3)
   { # @ARGV used in scalar context = number of args
	  print("This program uses apriori to mining the rules of GO terms in the uniprot database. \n");
	  print "Min support is set to 0.05, min confidence is set to 90!\n";
	  print "perl $0 addr_apriori addr_uniprot_dat dir_output\n";
	  print "For example:\n";
	  print "perl $0 /exports/store1/tool/apriori/src/apriori /exports/store1/tool/swiss_prot_2012/uniprot_sprot.dat ../result/R_1_mining_rules_from_uniprot\n";

	  exit(0);
   }
  

 my($addr_apri)=$ARGV[0];
 my($addr_uniprot)=$ARGV[1];
 my($dir_output)=$ARGV[2];
 -s $dir_output || system("mkdir $dir_output");
 -s $addr_apri || die "Cannot open the C program apriori, $addr_apri\n";
 -s $addr_uniprot || die "Cannot open the C program apriori, $addr_uniprot\n";

 my($IN,$OUT,$line,$path,$tmp);
 my(@tem_split,@tem_2);

 #######################################################################################################
 #
 #
 #                    1. parse uniprot database and get the GO terms
 #
 #
 #######################################################################################################
 $|=1;
 print "Start parsing the uniprot database ...";
 $path=$dir_output."/"."uniprot_parsed_GO_terms";             # store the parsed result, each line is the GO terms for one protein
 my($tag)=0;                                                  # to check whether the protein has GO term or not
 $OUT = new FileHandle ">$path";
 $IN=new FileHandle "$addr_uniprot";
 while(defined($line=<$IN>))
 {
	 chomp($line);
	 $line=~s/\s+$//;
	 @tem_split=split(/\s+/,$line);
	 if($line eq "//")
	 {
		 if($tag == 1)
		 {# this means we already print at least one GO term for this protein
			 print $OUT "\n";
		 }
		 $tag=0;
		 next;
	 }
	 if(($tem_split[0] eq "DR") && ($tem_split[1] eq "GO;"))
	 {
		 # we find one GO term for this protein
		 if($tag==0)
		 {# the first GO
			 @tem_2=split(/\;/,$tem_split[2]);

			 print $OUT $tem_2[0];
			 $tag=1;
		 }
		 else
		 {# we already find GO terms for this before
			 @tem_2=split(/\;/,$tem_split[2]);

			 print $OUT " ".$tem_2[0];
		 }
	 }
 }
 
 $IN->close();
 $OUT->close();

 #######################################################################################################
 #
 #
 #                    2. Use the tool to mining the rules 
 #
 #
 #######################################################################################################
 my($path_out)=$dir_output."/"."Mining_rules";
 my($return_val);
 $return_val = system("$addr_apri -tr -s0.05 -c90 $path $path_out");
 if($return_val!=0)
 {# we set min support as 0, support means the frequency of that rule exists in the whole database, it's OK if only exists few times, we need high confidence, confidence(x=>y) = support(X and Y)/Support(X). 
	 print "The program :$addr_apri -tr -s0.05 -c90 $path $path_out fails!\n";
 }
 else
 {
	 print "Great, done!\n";
 }

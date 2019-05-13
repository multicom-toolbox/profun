 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 sub make_prediction($);
   if (@ARGV != 3)
   { # @ARGV used in scalar context = number of args
	  print "This program makes protein function prediction from the sequence only! The GO is not filtered in this program. ALL GO terms is included\n";
	  print "perl $0 addr_new_Trained_sequence2GO addr_fasta_sequence addr_output_prediction\n";

	  print "For example:\n";
	  print "perl $0 ../result/new_Trained_sequence2GO ../data/CAFA2_all_fasta ../result/not_filter_GO_sequence_based_prediction_for_CAFA2\n";
          #print "perl $0 ../result/new_Trained_sequence2GO ../data/CAFA-2013-other/other-targets.fa ../result/not_filter_GO_sequence_based_prediction_for_other_species_CAFA2\n";

  	  exit(0);
   }

   my($s2g)=$ARGV[0];
   my($addr_seq)=$ARGV[1];
   my($output)=$ARGV[2];

   
 my($IN,$line,$OUT,$i);
 my(@tem_split,@tem_2,@tem_3);

 
 $|=1;
 print "reading sequence to GO information ...";
 my(%seq2GO)=();
 ##### read sequence to GO information #######
 $IN= new FileHandle "$s2g";
 while(defined($line=<$IN>))
 {
         chomp($line);
         @tem_split=split(/\s+/,$line); 
         if(@tem_split != 2)
	     {
			 print "Skip $line ..\n";
		 }
		 else
	     {
			 if(not exists $seq2GO{$tem_split[0]})
			 {
				 $seq2GO{$tem_split[0]}=$tem_split[1];
			 }
			 else
			 {
				 print "warning, duplicate for $line!\n";
			 }
		 }

 }
 $IN->close();
 print "Done!\n";
 
 print "Reading sequence and make prediction ...";
 ###### read sequence and make prediction #######
 my($OUT)=new FileHandle ">$output";

 my($sequence)="NULL";
 my($target_name)="NULL";
 my($predictions);
 $IN= new FileHandle "$addr_seq";
 while(defined($line=<$IN>))
 {
     chomp($line);
     if(substr($line,0,1) eq ">")
	 {
		 if($target_name ne "NULL")
		 {# we already get a sequence and target name, should make prediction on that 
			 if($sequence eq "NULL")
			 {
				 print "Error, this should not happen!\n";
				 exit(0);
			 }

             $predictions=make_prediction($sequence);         # GO_0.1|GO2_0.2 ...
			 if($predictions eq "NULL")
			 {# no prediction made 
				 print "no prediction for $sequence , next !\n";
				 next;
			 }
			 @tem_2=split(/\|/,$predictions);
			 for($i=0;$i<@tem_2;$i++)
			 {
				 @tem_3=split(/\_/,$tem_2[$i]);
				 print $OUT "Sequence_based\t$target_name\t$tem_3[0]\t$tem_3[1]\n";
			 }
		 }
		 @tem_split=split(/\s+/,$line);
		 @tem_split=split(/\>/,$tem_split[0]);
		 $target_name=$tem_split[1];
		 $sequence = "NULL";
	 }
     else
	 {
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split<1)
		 {# empty line
			 next;
		 }
		 if($sequence eq "NULL")
		 {
			 $sequence = $line;
		 }
		 else
		 {
			 $sequence.=$line;
		 }
	 }
     


 }
 ###### the last prediction before we finish ######
 $predictions=make_prediction($sequence);         # GO_0.1|GO2_0.2 ...
 if($predictions ne "NULL")
 {
  @tem_2=split(/\|/,$predictions);
  for($i=0;$i<@tem_2;$i++)
  {
	 @tem_3=split(/\_/,$tem_2[$i]);
	 print $OUT "Sequence_based\t$target_name\t$tem_3[0]\t$tem_3[1]\n";
  }
 }
 $OUT->close();

 sub make_prediction($)
 {# use %seq2GO to find the GO function of a sub sequence
	 my($seq)=@_;
     my($i1,$i2);
	 my($len)=length($seq);
	 my($sub_seq);
	 my(@temtem,@tt);
	 my(%hash)=();
     for($i1=0;$i1<$len-5;$i1++)
	 {
        $sub_seq=substr($seq,$i1,5);
		if(not exists $seq2GO{$sub_seq})
		{
			print "No $sub_seq, next!\n";
			next;
		}
		else
		{
            @temtem=split(/\|/,$seq2GO{$sub_seq});
			for($i2=0;$i2<@temtem;$i2++)
			{
				@tt=split(/\_/,$temtem[$i2]);

				if(not exists $hash{$tt[0]})
				{
					$hash{$tt[0]}=$tt[1];
				}
				else
				{
					$hash{$tt[0]}=1-(1-$hash{$tt[0]})*(1-$tt[1]);
				}
			}
		}
	 }
	 $sub_seq="NULL";
     foreach $i1 (keys %hash)
	 {
		 $hash{$i1} = sprintf("%.5f",$hash{$i1});
     }

	 foreach $i1 (sort {$hash{$b} cmp $hash{$a}} keys %hash)
	 {		 
		 
		 if($sub_seq eq "NULL")
		 {
			 $sub_seq = $i1."_".$hash{$i1};
		 }
		 else
		 {
			 $sub_seq .="|".$i1."_".$hash{$i1};
		 }
	 }
	 return $sub_seq;

 }

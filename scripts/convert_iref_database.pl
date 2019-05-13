 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 sub filter_file($$);
 use POSIX qw/ceil/; #this is for function ceil()
  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
	  print("This program file filters the irefindex database, only keep the ppi networks!\n");
	  print "perl $0 dir_input dir_output\n";
          print "For example:\n";
          print "perl $PROGRAM_NAME ../data/ireindex_orig_database ../data/ireindex_converted_database\n";
      exit(0);
	}
 my($dir_input)=$ARGV[0];
 my($dir_output)=$ARGV[1];
 
 -s $dir_input || die "Cannot find input file $dir_input!\n";
 my($IN,$line,$file);
 my(@tem_split,@files);
 my($path_input,$path_output);
 -s $dir_output || system("mkdir $dir_output");
 opendir(DIR,$dir_input);
 @files=readdir(DIR);
 foreach $file (@files)
 { 
    if($file eq '.' || $file eq '..')
    {
       next; 
    }
    print "processing $file ...\n";
    $path_input=$dir_input."/".$file;
    $path_output=$dir_output."/".$file;
    filter_file($path_input,$path_output); 
 }
 
 sub filter_file($$)
 {
    my($input,$output)=@_;
    my($line,$ID1,$ID2);

    ########### we don't want duplicate edeges ############
    my(%hash_node)=();
    my($node);
    #######################################################

    my(@tem_split,@tem2);
    my($IN)=new FileHandle "$input";
    my($OUT) = new FileHandle ">$output";
    defined($IN) || die "Cannot open $input!\n";
    while(defined($line=<$IN>))
    {
         chomp($line);
         $line=~s/\s+$//;
         @tem_split=split(/\s+/,$line);
         if(@tem_split<1)
         {
                 next;
         }
         @tem2 = split(/\:/,$tem_split[0]);
         if($tem2[0] eq "uniprotkb")
         {
              $ID1=$tem2[1];
              @tem2=split(/\:/,$tem_split[1]);
              if($tem2[0] eq "uniprotkb")
              {
                  $ID2=$tem2[1];
                  if($ID1 eq $ID2)
                  { # we skip the self contact
                     #print $line."\n";
                     next;
                  }
#print $line."\n";
                  ######## check whether this edge exists ##########
                  $node=$ID1."_".$ID2;
                  if(exists $hash_node{$node})
                  {
                      next;
                  }
                  else
                  {
                      $node=$ID2."_".$ID1;
                      if(exists $hash_node{$node}) {next;}
                  }
                  $hash_node{$node}=1;
                  ###### we don't keep the ID with - inside, for exampe A0JLT2-2 is used as A0JLT2 #####
                  @tem2=split(/\-/,$ID1);
                  if(@tem2>1)
                  {
                     $ID1=$tem2[0];
                  }
                  @tem2=split(/\-/,$ID2);
                  if(@tem2>1)
                  {
                     $ID2=$tem2[0];
                  }
                  ######################################################################################
                  print $OUT $ID1."\t".$ID2."\n";
              }
         } 
         
    }
    $IN->close();
    $OUT->close();
 }

 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 sub filter_file($$);
 sub parse_separate_network($$$);       # parse and get separate network for one specific species
 use POSIX qw/ceil/; #this is for function ceil()
  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
          print "Updating at 10/5/2013, only keep the protein protein interaction in the same sepcies, delete the protein interacting with other proteins \n";
      print "Script just for CAFA2, I want to find all CAFA2 species network, I am not sure whether I can get that!\n";

	  print("This program file filters the irefindex database of All networks, get the networks with the taxid\n");
	  print "perl $0 addr_input(it's the All network of irefindex) dir_output\n";
          print "For example:\n";
          print "perl $PROGRAM_NAME ../download_database/All.mitab.06062013.txt ../data/converted_network_with_taxid\n";
      exit(0);
	}
 my($addr_input)=$ARGV[0];
 my($dir_output)=$ARGV[1];
 

 my($IN,$line,$file,$OUT);
 my(@tem_split,@files);
 my($i,$path_input,$i,$path_output);
 -s $dir_output || system("mkdir $dir_output");
 
 $path_input=$addr_input;
 @tem_split=split(/\//,$addr_input);
 
 my(@all_ids)=();
 $path_output=$dir_output."/".$tem_split[@tem_split-1];

 
 @all_ids=filter_file($path_input,$path_output);  
 
 print "All ids are as follows:\n";
 print "@all_ids"."\n";
 
 
 @all_ids=();
 $all_ids[0] = "186497";
 $all_ids[1] = "243232";
 $all_ids[2] = "273057";
 $all_ids[3] = "309800";
 $all_ids[4] = "436308";
 $all_ids[5] = "453591";
 $all_ids[6] = "478009";
 $all_ids[7] = "160488";
 $all_ids[8] = "170187";
 $all_ids[9] = "208964";

 $all_ids[10] = "223283";
 $all_ids[11] = "224308";
 $all_ids[12] = "243273";
 $all_ids[13] = "321314";
 $all_ids[14] = "83333";
 $all_ids[15] = "85962";
 $all_ids[16] = "99287";
 $all_ids[17] = "10090";
 $all_ids[18] = "10116";
 $all_ids[19] = "284812";
 $all_ids[20] = "3702";
 $all_ids[21] = "44689";
 $all_ids[22] = "559292";
 $all_ids[23] = "7227";
 $all_ids[24] = "7955";
 $all_ids[25] = "8355";
 $all_ids[26] = "9606"; 
 $all_ids[27] = "other";

 $path_input=$path_output;
 
 my(%hash_all_ids)=();
 for($i=0;$i<@all_ids;$i++)
 {
	 $hash_all_ids{$all_ids[$i]}=$i;
 }
 ##### first delete existing result #######
 for($i=0;$i<@all_ids;$i++)
 {
	 $path_output=$dir_output."/".$all_ids[$i].".network";
	 if(-s $path_output)
	 {
		 system("rm $path_output");
	 }
 }
 
 $IN=new FileHandle "$path_input";
 defined($IN) || die "Cannot open $path_input!\n";
 while(defined($line=<$IN>))
 {
      chomp($line);
      $line=~s/\s+$//;
      @tem_split=split(/\s+/,$line);
	  if($tem_split[2] eq $tem_split[3])
	  {
			 if(exists $hash_all_ids{$tem_split[2]})
			 {
			      $path_output=$dir_output."/".$all_ids[$hash_all_ids{$tem_split[2]}].".network";
				  $OUT = new FileHandle ">>$path_output";
				  print $OUT $line."\n";
			      $OUT->close();
			 }
			 else
		     {
			      $path_output=$dir_output."/".$all_ids[27].".network";
				  $OUT = new FileHandle ">>$path_output";
				  print $OUT $line."\n";
			      $OUT->close();
			 }
	  }
	  else
	  {# different species
####### we don't consider interaction between species
=pod
			  if(exists $hash_all_ids{$tem_split[2]})
			  {
                  $path_output=$dir_output."/".$all_ids[$hash_all_ids{$tem_split[2]}].".network";
				  $OUT = new FileHandle ">>$path_output";
				  print $OUT $line."\n";
				  $OUT->close();
			  }
			  if(exists $hash_all_ids{$tem_split[3]})
			  {
                  $path_output=$dir_output."/".$all_ids[$hash_all_ids{$tem_split[3]}].".network";
				  $OUT = new FileHandle ">>$path_output";
				  print $OUT $line."\n";
				  $OUT->close();
			  }
			  if((not exists $hash_all_ids{$tem_split[2]}) && (not exists $hash_all_ids{$tem_split[3]}))
		      {
			      $path_output=$dir_output."/".$all_ids[8].".network";
				  $OUT = new FileHandle ">>$path_output";
				  print $OUT $line."\n";
			      $OUT->close();
			  }
=cut
	  }


 }

 $IN->close();

# for($i=0;$i<@all_ids;$i++)
# {# filter each separate species, and get specific network !
#    $path_output=$dir_output."/".$all_ids[$i].".network";
#	parse_separate_network($path_input,$path_output,$all_ids[$i]);
# }
 
 #@all_ids=

=pod
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
=cut 
 sub parse_separate_network($$$)
 {
	 my($input,$output,$id)=@_;
	 my($IN,$OUT,$line);
	 my(@tem_split);
    my($IN)=new FileHandle "$input";
    my($OUT) = new FileHandle ">$output";
    defined($IN) || die "Cannot open $input!\n";
    while(defined($line=<$IN>))
    {
         chomp($line);
         $line=~s/\s+$//;
         @tem_split=split(/\s+/,$line);
		 if($tem_split[2] eq $id || $tem_split[3] eq $id)
		 {
			 print $OUT $line."\n";
		 }
	}
	$IN->close();
	$OUT->close();
 }
 sub filter_file($$)
 {
    my($input,$output)=@_;
    my($line,$ID1,$ID2);

    ########### we don't want duplicate edeges ############
    my(%hash_node)=();
    my($node);
    #######################################################
    
    my(%hash_ids)=();
	my(@ids)=();
	my($i,$key,$id1,$id2);


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

				  
				  ####### add the new taxids #####

				  @tem_split=split(/\t/,$line);

				  if($tem_split[9] eq "-" || $tem_split[10] eq "-")
				  {
					  print "$line, taxid is $tem_split[9] and $tem_split[10]\n";
					  #exit(0);
				  }

				  @tem2=split(/\(/,$tem_split[9]);
				  if(not exists $hash_ids{$tem2[0]})
				  {
					  $hash_ids{$tem2[0]}=1;
				  }
				  $key=$tem2[0];
				  @tem2=split(/\:/,$tem2[0]);
				  
				  $id1=$tem2[1];
				  @tem2=split(/\(/,$tem_split[10]);
				  if(not exists $hash_ids{$tem2[0]})
				  {
					  $hash_ids{$tem2[0]}=1;
				  }
				  #if($key ne $tem2[0])
				  #{
					  #print "at line  connection between two species: $key and $tem2[0]\n";
				  #}
				  @tem2=split(/\:/,$tem2[0]);
				  $id2=$tem2[1];
				  ################################
				  ###### we don't use the ID with - inside ########
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
                  print $OUT $ID1."\t".$ID2."\t".$id1."\t".$id2."\n";
              }
         } 
         
    }
    $IN->close();
    $OUT->close();
    $i=0;
	foreach $key (keys %hash_ids) 
    {
		$ids[$i]=$key;
		$i++;
	}
	return @ids;
 }

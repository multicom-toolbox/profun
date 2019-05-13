$num = @ARGV;
if($num != 3)
{
	die "The number of parameter is not correct!\n";
}

$input_dir = $ARGV[0];
$output_dir = $ARGV[1];
$scripts =  $ARGV[2];

opendir(DIR,"$input_dir") || die "Failed to open file $input_dir\n";
@files=readdir(DIR);

my($total_number)=scalar(@files);
foreach $file (@files) 
{
      if($file eq "." || $file eq ".." || index($file,'network')<0)
      {
           next;
      }

      if($file eq '559292.network' or $file eq '10090.network' or $file eq '9606.network' or $file eq '10116.network' or $file eq '7227.network' or $file eq '284812.network' or $file eq '83333.network' )
      {
	print "10090.network and 559292.network is running by Jie, pass \n";
	next;
      }

      $count_number++;
      @temp = split(/\./,$file);
      $name = $temp[0];
      $networkfile = $input_dir.'/'.$file;
      print "processing $networkfile\n";
      $outputfile = $output_dir.'/'.$name.'.go_frequency';
      if(-e $outputfile)
      {
          print "$outputfile has been generated, pass \n";
	  next;
      }
      print("perl $scripts/A_3_cal_go_pairs_frequency.pl  $networkfile $outputfile");
      $status = system("perl $scripts/A_3_cal_go_pairs_frequency.pl  $networkfile $outputfile");

      if($status)
      {
	die "Failed to run <perl $scripts/A_3_cal_go_pairs_frequency.pl  $networkfile $outputfile>\n"
      }

}



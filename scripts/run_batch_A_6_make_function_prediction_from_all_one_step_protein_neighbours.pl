$num = @ARGV;
if($num != 6)
{
        die "The number of parameter is not correct!\n";
}

$go_frequency_input_dir = $ARGV[0];
$real_go_function_input_dir =  $ARGV[1];
$one_step_protein_neighbour_input_dir =  $ARGV[2];
$real_go_function_input_dir2 =  $ARGV[3];
$output_dir = $ARGV[4];
$scripts =  $ARGV[5];

opendir(DIR,"$go_frequency_input_dir") || die "Failed to open file $go_frequency_input_dir\n";
@files=readdir(DIR);

my($total_number)=scalar(@files);
foreach $file (@files)
{
      if($file eq "." || $file eq ".." || index($file,'go_frequency')<0)
      {
           next;
      }

      $count_number++;
      @temp = split(/\./,$file);
      $name = $temp[0];
      $go_frequency_file = $go_frequency_input_dir.'/'.$file;
      $real_go_function_file = $real_go_function_input_dir.'/'.$name.'.mitab.go_functions';
      $one_step_protein_neighbour_file = $one_step_protein_neighbour_input_dir.'/'.$name.'.one_step_protein_neighbours';
      $real_go_function_file2 = $real_go_function_input_dir2.'/'.$name.'.mitab.go_functions';
      print "processing $go_frequency_file\n";
      $outputfile = $output_dir.'/'.$name.'.function_prediction';
      if(-e $outputfile)
      {
          print "$outputfile has been generated, pass \n";
          next;
      }
      print("perl $scripts/A_6_make_function_prediction_from_all_one_step_protein_neighbours.pl $go_frequency_file $real_go_function_file $one_step_protein_neighbour_file $real_go_function_file2      $outputfile");
      $status = system("perl $scripts/A_6_make_function_prediction_from_all_one_step_protein_neighbours.pl $go_frequency_file $real_go_function_file $one_step_protein_neighbour_file $real_go_function_file2 $outputfile");

      if($status)
      {
        die "Failed to run <perl $scripts/A_6_make_function_prediction_from_all_one_step_protein_neighbours.pl $go_frequency_file $real_go_function_file $one_step_protein_neighbour_file $real_go_function_file2  $outputfile>\n"
      }

}

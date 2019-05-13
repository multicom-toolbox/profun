import java.util.Stack;  
import java.io.*;
import java.util.ArrayList;


public class TestSequenceAll {  
  public static void main(String[] args) {  
    TestSequenceAll t = new TestSequenceAll();

    //Object[] arr = { "GO:1", "GO:2", "GO:3" , "GO:4" , "GO:5" , "GO:6", "GO:7", "GO:8", "GO:9", "GO:10", "GO:11", "GO:12"}; 
	//ArrayList arr = new ArrayList();
	Object[] arr;
	arr = new Object[200];

	//System.out.println(args.length);
	if(args.length<2)
	{
		System.out.println("Error!\nThe input parameter should be path_of_input_data and total_number_of_combination_bits ");
		System.exit(0);
	}
	int index4arr=0;
	int mallocLen=0;

    //////////////  Read the input args[0], and store it in arr /////////////////
        File file = new File(args[0]);
        BufferedReader reader = null;
        try {
            
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;

            while ((tempString = reader.readLine()) != null) {
				if(tempString.equals(""))
				{
					continue;
				}
                arr[index4arr]=tempString;
				index4arr++;
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
        }
	/////////////////////////////////////////////////////////////////////////////
	int i;

    int totalOut = index4arr;
	if(index4arr > Integer.parseInt(args[1]))
	{
		totalOut=Integer.parseInt(args[1]);
	}
	else
	{
		totalOut = index4arr;
	}
//    System.out.println("The totalOut is set to " + totalOut);

    for (int num = 1; num <= totalOut; num++) {  
		
        t.getSequence(arr, 0, num, index4arr);  
    }  
  }  

  private Stack<Object> stack = new Stack<Object>();  
  public void getSequence(Object[] arr, int begin, int num, int index4arr) {  
    if (num == 0) {  
      System.out.println(stack);
    } else {  

      for (int i = begin; i < index4arr; i++) {  

        stack.push(arr[i]);  

        swap(arr, begin, i);  

        getSequence(arr, begin + 1, num - 1, index4arr);  

        swap(arr, begin, i);  

        stack.pop();  
      }  
    }  
  }  

  public static void swap(Object[] arr, int from, int to) {  
    if (from == to) {  
      return;  
    }  
    Object tmp = arr[from];  
    arr[from] = arr[to];  
    arr[to] = tmp;  
  }  
}

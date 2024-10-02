template<typename t>
 t getMax(t *arr, int n) 
{
    t max = arr[0];
   for (int i = 0; i < n; i++) 
   {
    if (max < arr[i]) 
	{
     max = arr[i];
	}
   }
     return max;
 }

 template<typename t>
 t getMin(t *arr, int n)
 {
	 t max = arr[0];
	 for (int i = 0; i < n; i++)
	 {
		 if (max > arr[i])
		 {
			 max = arr[i];
		 }
	 }
	 return max;
 }




 //template<typename t>
 //t Pick(t *arr, int start, int end)
 //{
	// int num = end - start + 1;
	// t aa = new t[num]();
	// for (int i = 0; i < end - start + 1; i++)
	// {
	//	 aa[i] = arr[start+i-1];
	// }
	// return aa;
 //}

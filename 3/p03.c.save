#include <stdio.h>
#include <math.h>
#include <string.h>

int main()
{

char seq[5000];
int len;
float noCGs;
int x;




while (scanf("%s", seq)==1)
{
	
	len = strlen( seq );

	for ( x =0 ; x < len; x++)
	 	{
		if(seq[x] == 'C' || seq[x] =='c'|| seq[x] == 'G' || seq[x] =='g')
			{
			noCGs ++;
			}
		}


	printf("The sequence is %s, and its length is %d bp.", seq, len);
	printf("The percentage of GCs is %f %%. \n", noCGs*100/len);
}	


}


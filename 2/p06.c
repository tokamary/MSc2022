#include <stdio.h>
#include <math.h>

int main()
{

char k;
float	total;
float seqlen;
float y;

total = 0;
seqlen = 0;
y = 0.00;


	while( scanf("%c", &k )==1)
	{
	seqlen ++;

	if ( k == 'C'|| k=='c'||k=='G'||k=='g')
		{
		total ++;
		}
	}
	
	y=total*100/seqlen;

	printf("The persentage of Cs and Gs is %f. %% \n",y);

}

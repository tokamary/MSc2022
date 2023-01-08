#include <stdio.h>
#include <math.h>

int main()
{

char k;
int	total;
total = 0;

	while( scanf("%c", &k )==1)
	{
	if ( k == 'C'|| k=='c')
		{
		total ++;
		}
		
	}
	
	printf("The total number of Cytosine is %d.\n", total);

}

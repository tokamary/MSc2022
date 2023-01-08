#include <stdio.h>
#include <math.h>

int main()
{

float	x;
float	y;
x = 0.0;

while (x <= 100.0)
	{
		y = sqrt(x);
		printf("The square root of %f is %f. \n", x, y);
	 	x= x + 3;

	}



}


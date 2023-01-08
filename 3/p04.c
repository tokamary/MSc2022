#include <stdio.h>
#include <math.h>
#include <string.h>

int main()
{

char seq[5000];
int len;
int x;

while (scanf("%s", seq)==1)
{
	len=strlen(seq);

	for(x=0 ; x < len ; x++)
	{
	if(seq[x]=='D' || seq[x]=='E' || seq[x]=='H' ||seq[x]=='K' ||seq[x]=='N' ||seq[x]=='Q'|| seq[x]=='R')
		{
		printf(" ");
		}
	if(seq[x]=='S' || seq[x]=='T' || seq[x]=='G')
		{
		printf(".");
		}
	if(seq[x]=='A' || seq[x]=='C' || seq[x]=='M' ||seq[x]=='P')
		{
		printf(":");
		}
	if(seq[x]=='F' || seq[x]=='I' || seq[x]=='L' ||seq[x]=='V' ||seq[x]=='W' ||seq[x]=='Y')
		{
		printf("*");
		}
	}
}




}


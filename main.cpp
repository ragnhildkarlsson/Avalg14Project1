#include <iostream>
#include <stdio.h>
#include <gmp.h>

using namespace std;




void printBigInt(){
	mpz_t integ;
	mpz_init (integ);
	mpz_set_str(integ,"2849585738957348573897639763487638376",10);
  	printf ("integ");
  	mpz_out_str(stdout,10,integ);
  	printf ("\n");
  	mpz_clear(integ);
}


int main()
{

	printBigInt();
	return (0);
}






#include <iostream>
#include <stdio.h>
#include <gmpxx.h>


using namespace std;




void printBigInt(){
	mpz_class a;

	a = "1284497587348578763896789346950";
	cout <<  "a is " << a << "\n";

}


int main()
{

	printBigInt();
	return (0);
}




	// mpz_init (integ);
	// mpz_set_str(integ,"2849585738957348573897639763487638376",10);
 //  	printf ("integ");
 //  	mpz_out_str(stdout,10,integ);
 //  	printf ("\n");
  	// mpz_clear(integ);


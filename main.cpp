#include <iostream>
#include <stdio.h>
#include <gmpxx.h>


using namespace std;




void printBigInt(){
	mpz_class a;

	a = "1284497587348578763896789346950";
	cout <<  "a is " << a << "\n";

}


int isPrime(mpz_class n){
	
	//
	mpz_class rest, t, s, a, test;
	gmp_randclass rand (gmp_randinit_default);

	if(n==2){
		return 1;
	}
	//check if n is even
	rest = n % 2;
	if(rest == 0){
		return 0;
	}

	//calculate t as n-1 = 2^s * t
	t = n-1;
	s=0;
	while(t%2 == 0){
		t= t >> 1;
		cout << "t is "<< t<< "\n";
		s=s+1;
	}
	cout << "s is "<< s<< "\n";

	//a is a random number (1<=a<=N-1)
	a = rand.get_z_range(n-1) +1;

	//we compute a^t mod 

	

	//print binary representation of 67
	test = 67;
	cout << "bin-test is "<< test.get_str(2) << "\n" ;

	//we comput a 
	//find u0 such as u0 congruent were a is a random number (1<=a<=N-1)


	return 1;
}

mpz_class modExponential(mpz_class base, mpz_class exp, mpz_class mod){
		mpz_class res;
		mpz_t rop;
		mpz_init(rop);
		char c;
		string binExponent = exp.get_str(2);
		mpz_set(rop, base.get_mpz_t());

		for(unsigned i=1; i< binExponent.length();i++){			
			c = binExponent.at(i);
			
			if(c=='1'){

				mpz_pow_ui(rop, rop, 2);
				res = mpz_class(rop);
				res = (res*2) % mod;
				mpz_set(rop, res.get_mpz_t());
			}
			if(c=='0'){
				mpz_pow_ui(rop, rop, 2);
				res = mpz_class(rop);
				res = (res) % mod;
				mpz_set(rop, res.get_mpz_t());

			}

		}
		return res;
}


int main()
{
	mpz_class base, exp, mod;
	base= 2;
	exp = 154;
	mod =155;
		
	cout<< modExponential(base, exp,mod) <<"\n";
	return (0);
}




	// mpz_init (integ);
	// mpz_set_str(integ,"2849585738957348573897639763487638376",10);
 //  	printf ("integ");
 //  	mpz_out_str(stdout,10,integ);
 //  	printf ("\n");
  	// mpz_clear(integ);


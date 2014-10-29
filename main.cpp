#include <iostream>
#include <stdio.h>
#include <gmpxx.h>


using namespace std;




void printBigInt(){
	mpz_class a;

	a = "1284497587348578763896789346950";
	cout <<  "a is " << a << "\n";

}

mpz_class pow(mpz_class base, unsigned long exp){
		mpz_class res;
		mpz_t temp;
		mpz_init(temp);
		
		mpz_set(temp, base.get_mpz_t());
		mpz_pow_ui(temp, temp, exp);
		res = mpz_class(temp);
		mpz_clear(temp);
		return res;
		
}

mpz_class modExponential(mpz_class base, mpz_class exp, mpz_class mod){
		mpz_class res;
		res = base;
		char c;
		string binExponent = exp.get_str(2);
		for(unsigned i=1; i< binExponent.length();i++){			
			c = binExponent.at(i);

			if(c=='1'){
				res = pow(res, 2);
				res = (res*base) % mod;
				
			}
			if(c=='0'){
				res = pow(res,2);
				res = (res) % mod;
			}
			cout << "res is "<< res<< "\n";
		}
		return res;
}


int isPrime(mpz_class n, int rep){
	
	mpz_class rest, t, s, a, u;
	gmp_randclass rand (gmp_randinit_default);
	bool isPrime;
	isPrime = false;	
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
		//cout << "t is "<< t<< "\n";
		s=s+1;
	}
	//cout << "s is "<< s<< "\n";

	for(int i =0; i< rep;++i){
		//a is a random number (1<=a<=N-1)
		a = rand.get_z_range(n-1) +1;
		//cout << "\n Rand value a is " << a;
		//cout << "\n value of t " << t << " value of s " << s;
		//we compute a^t mod 
		u= modExponential(a,t,n);
		//cout << "\n u0 is " << u;
		if(u==1 || u == n-1){
			isPrime = true;
		}
		else{
			for(int j =1;j<=s;++j ){
				u = pow(u,2);
				//cout << "\n u" << j << " is " << u;
				if(u==n-1){
					isPrime =true;
					break;
				}
			}

		}
		if(!isPrime){
			return 0;
		}
		isPrime=false;
	}

	return 1;
}



int main()
{
	mpz_class t1, h1,h2,h3;
	t1 = "169743212304909";
	h1 = 5;
	h2 = 3;
	h3 = 7;
	
	cout << "\n t1 result "<< isPrime(t1, 10)<< "\n";

		
	return (0);
}




	// mpz_init (integ);
	// mpz_set_str(integ,"2849585738957348573897639763487638376",10);
 //  	printf ("integ");
 //  	mpz_out_str(stdout,10,integ);
 //  	printf ("\n");
  	// mpz_clear(integ);


#include <iostream>
#include <stdio.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <limits>



using namespace std;


static const  int primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};
static const int primeLength = 26;

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
mpz_class gcd(mpz_class n1, mpz_class n2){
	mpz_class t, u,v ;
	if(n1>n2){
		u = n1;
		v= n2;
	}
	else{
		u = n2;
		v= n1;
	}
	while(v>0){
		t= u;
		u=v;
		v= t%v;
	}

	return u <0 ? -u : u;
}


mpz_class * genTestData(int size, int j){
	
	mpz_class sufix,ten, startval;
	ten = 10;

	sufix = pow(ten, 60+j);
	static mpz_class data[200];
	string micke = "9207064750";
	string ragnhild = "8502142782";
	startval = micke;
	startval = startval*sufix;
	int index = 0;
	for(int i=1; i<= size;++i){
		data[index] = startval+i;
		++index;
	}
	startval = ragnhild;
	startval = startval*sufix;	
	for(int i=1; i<= size;++i){
		data[index] = startval+i;
		++index;
	}
	return data;


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
			//	cout << "res is "<< res<< "\n";
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
				u = u % n;
		//		cout << "\n u" << j << " is " << u;
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

mpz_class pollard(mpz_class seed, mpz_class add, mpz_class N){
	// check small primes

	for(int i = 0; i < primeLength; ++i){
		if((N % primes[i]) == 0){
			return primes[i];
		}
	}

	mpz_class a,b, p, product, medium, maxint;
	maxint =  std::numeric_limits<int>::max();
	product = 1;
	a = seed;
	b = seed;

	int maxIter;
	int maxInt;
	maxIter = 1000;

	// //Large N
	maxInt = std::numeric_limits<int>::max();
	// maxIter = maxInt;

	// //Medium N
	// medium =10;
	// medium = pow(medium,20);

	// if(N<medium){
	// 	maxIter = 100;	
	// }
	
	//small N
	if(N<maxInt){
		maxIter = 1;
	}
	
	do{
		for(int i =0; i<maxIter;i++){
		a=(pow(a,2) - add) % N;
		b=(pow(b,2) -add) % N;
		b =(pow(b,2) -add) % N;
		b =(pow(b,2) -add) % N;
		//cerr << "\n a is " << a << " b is " << b;	

		if(a == b) {
			break;

		}
		product = product*abs(a-b);
		
		}
		p = gcd(product,N);
		if(p>1){
			return p;
		}

	}while(a!=b);

	p = "-1";
	return p;

}


void factorize(mpz_class N, int trials){
	cout << "\n Trying to factorize " << N;
	
	std::clock_t start;
    double duration;

    start = std::clock();

    cout << start;
	vector<mpz_class> factors;	
	gmp_randclass rand (gmp_randinit_default);
	
	mpz_class seed, add, factor, temp;
	bool primFound = false;
	while(N!=1){
		factor = N;
		//cerr << "\nFactor in the outer loop is now " << factor;		
		if(isPrime(factor,10)){
			factors.push_back(factor);
			cerr << "\n Added factor " << factor;
			break;
		}

		while(!primFound){
			cerr << "\n Entering primFound loop with factor " << factor;
			for(int i = 0; i < trials; ++i){
				seed = rand.get_z_range(factor-1) +1;
				add = 1; //TODO Hard CODED
				cerr << "\n Seed is " << seed << " add is " << add;
				temp = pollard(seed, add, factor);
					if(temp > 0){
						factor = temp;
						//cerr << "\nIs the factor " << factor << " considered prime? " << isPrime(factor,10);
						if(factor > 0 && isPrime(factor,10)){
							factors.push_back(factor);
							cerr << "\n Added factor " << factor;
							primFound = true;
							break;
						}
						// Compsite number found, try to factor it into a prime
						if(factor > 0) {
							break;
						}
					}
				}
				// Tried trial times and still failed?
				if(temp < 0){
					cout << "\n Factorization failed!";
					return;
				}

			}
			N = N / factor;
			primFound = false;

		}
	cout << "\n The number of prime factors is " << factors.size(); 
	//print factors
	for (unsigned i = 0; i < factors.size() ; ++i)
	{
		cout <<"\n " << factors.at(i);
	}
	cout << "\n";
	duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;
    cout<<"Calculation time: "<< duration <<'\n';	
	}



int main()
{

	int trials;
	trials = 2;

	 mpz_class h1,h2, h3,t1;
	 h1 = 561;
	 h2 = 1;
	 h3 = 1729;
	 //1729

	 mpz_class* data = genTestData(100,20);
	 //data 4 bra exempel på pollard 
	 t1 = data[1];
	 cout << "\nt1 is " << t1;
	 //cout << "\npollard says " << pollard(h1,h2,h3);


	 //	cout << "is it red alert here?";
	 factorize(t1, trials);
	 //factorize(h2, trials);	

	// cout << "test 1 " <<gcd(h1,h2);
	// cout << "test 2 " <<gcd(h2,h1);
	// cout << "test 1 " <<pollard(h1,h2,h3);	
	// cout << "\n t1 result "<< isPrime(t1, 10)<< "\n";
	// 
	// for(int i=0; i<100*2;++i){
	// 	cout << "\n " << data[i];
	// }

		
	return (0);
}




	// mpz_init (integ);
	// mpz_set_str(integ,"2849585738957348573897639763487638376",10);
 //  	printf ("integ");
 //  	mpz_out_str(stdout,10,integ);
 //  	printf ("\n");
  	// mpz_clear(integ);


#include <iostream>
#include <stdio.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <limits>
#include <math.h>


using namespace std;


static const  int smallPrimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};
static const int primeLength = 26;
static const int primeSafety = 10;
static const int maxB = 5000000;
static const double e = 2.718281;



mpz_class pow(mpz_class base, unsigned long exp){
		mpz_class res;
		mpz_t temp;

		if(exp ==0){
			res =1;
			return res;
		}
		if( exp ==1){
			return base;
		}

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


mpz_class * genTestData(int j){
	
	mpz_class sufix,ten, startval;
	ten = 10;

	sufix = pow(ten, 60+j);
	static mpz_class data[200];
	string micke = "9207064750";
	string ragnhild = "8502142782";
	startval = micke;
	startval = startval*sufix;
	int index = 0;
	for(int i=1; i<= 100;++i){
		data[index] = startval+i;
		++index;
	}
	startval = ragnhild;
	startval = startval*sufix;	
	for(int i=1; i<= 100;++i){
		data[index] = startval+i;
		++index;
	}
	return data;


}


mpz_class modExponential(mpz_class base, mpz_class exp, mpz_class mod){


		mpz_class res;
		res = base;
		char c;
		if(exp ==0){
			res =1;
			return res;
		}

		if(exp ==1){
			res = res % mod;
			return res;	
		}
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
				u = modExponential(u,2,n);
				
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

	//std::clock_t pClock;
    //double duration;
    //pClock = std::clock();

	for(int i = 0; i < primeLength; ++i){
		if((N % smallPrimes[i]) == 0){
			return smallPrimes[i];
		}
	}

	mpz_class a,b, p, product, medium, maxint;
	double doubleN ,lowerBound, upperBound, bound;
	doubleN = N.get_d();
	lowerBound = log(doubleN);
	upperBound = sqrt(sqrt(doubleN));
	bound = (upperBound+lowerBound)/2;
	bound = lowerBound;
		cerr <<"\n upperBound" <<upperBound;
		cerr << "\n using bound "<< bound;
	
	product = 1;
	a = seed;
	b = seed;

	
	do{
		for(double i =0; i<bound;i++){
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
		//duration = (std::clock() - pClock ) / (double) CLOCKS_PER_SEC;
	}while(a!=b);


	
    //cout<<"Calculation time: "<< duration <<'\n';	



	p = "-1";
	return p;

}
 std::vector<mpz_class> genPrimeBase(){
	mpz_class i;
	i =2;
	vector<mpz_class> primes;
	primes.push_back(i);
	i=3;
	while(i < maxB){
		if(isPrime(i, primeSafety)){

			primes.push_back(i);
		}

		i = i+2;
	}
 return primes;

}
int legendreSymbol(mpz_class a, mpz_class p){
	mpz_class temp;
	if(a % p ==0){
		return 0;
	}
	unsigned long exponent = (p.get_ui() -1)/2;
	temp = modExponential(a, exponent, p);
	if(temp==1){
		return 1;
	}
	return -1;

}

mpz_class shanksTonelli(mpz_class n,mpz_class p){
	mpz_class q ,s,m, z,t,b,c,two,temp, res;
	if(p==2){
		res = 2;
		return res;
	}
	if(p%4 ==3){
		res = modExponential(n, (p+1)/4, p);
		return res;
	}

	//calculate t as n-1 = 2^s * t
	q = p-1;
	s=0;
	while(q%2 == 0){
		q= q >> 1;
		//cout << "t is "<< t<< "\n";
		s=s+1;
	}
	//Find z = a non-quadratic residue in p
	z=2;
	while(legendreSymbol(z,p)!=-1){
		n=n+1;
	}
	c = modExponential(z,q,p);
	res = modExponential(n, (q+1)/2, p);
	t = modExponential(n,q,p);
	m=s;

	while(true){
		cerr <<"\nInside infinity";
		cerr << "\n t: " <<t <<" c: "<< c << " res " << res <<" m: "<< m<<" z: "<< z <<" q: "<< q;
		if(t%p ==1){
			return res;
		}
		int i;
		temp=t;
		for (i = 1; i < m; ++i)
		{
			temp = modExponential(temp,2,p);
			if(temp%p ==1){
				break;
			}
		}
		cerr <<"\n after for loop is i: "<<i <<"is m "<<m <<" is t "<< t;

		two =2;

		unsigned long expo = m.get_ui();
		expo = expo-1-i; 
		cerr <<"\n expo is " << expo;

		b= modExponential(c, pow(two,expo) , p);
		res= (res*b)  % p;

		t = (t*b*b) % p;
		c= (b*b)%p;
		m = i;
		cerr << "\n b: " <<b <<" c: "<< c << " res " << res <<" m: "<< m<<" t: "<< t;
		
	}
}

int quadraticSieve(mpz_class N, vector<mpz_class> & primes){
	
	mpz_class temp;
	cout << "\n QS trying to factorize \n" << N <<"\n";
	vector<mpz_class> factorBase;
	vector<mpz_class> Qx;


	//Create factorbase
	double B, doubleN, bExp;
	doubleN = N.get_d();
	bExp = sqrt(log(doubleN)*log(log(doubleN)))*0.5;
	B = 3*pow(e, bExp);
	temp = 2;
	factorBase.push_back(temp);
	
	B = 100;
	for (unsigned i = 1; i < primes.size(); ++i)
	{	
		if(primes[i]>B){
			break;
		}
		unsigned long exponent = (primes.at(i).get_ui() -1)/2;
		temp = modExponential(N, exponent, primes.at(i));
		if(temp == 1){
			factorBase.push_back(primes.at(i));
		} 
	}
	//	Factorbase done

	//Create plynominals
	// log p, t, -t (//implement shanks-tonelli)



	return 0;

}


int factorize(mpz_class N, int trials){
	cout << "\n Trying to factorize \n" << N <<"\n";
	
	std::clock_t start;
    double duration;

    start = std::clock();

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
			cerr << "\n Entering primFound loop with factor: \n " << factor << "\n Trial nr: ";
			for(int i = 0; i < trials; ++i){
				seed = rand.get_z_range(factor-1) +1;
				add = 1; //TODO Hard CODED
				cerr << i << " ";
				//cerr << "\n Seed is \n " << seed;
				temp = pollard(seed, add, factor);
					if(temp > 0){
						factor = temp;
						//cerr << "\nIs the factor " << factor << " considered prime? " << isPrime(factor,10);
						if(factor > 0 && isPrime(factor,primeSafety)){
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
					return 0;
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
    cout<<"\n Calculation time: "<< duration <<'\n';	
	return 1;

	}



int main()
{	

	// cerr << "\nGenerating primes...";
	// vector<mpz_class> primes;
	
	// primes = genPrimeBase();
	// cerr <<"\nsize of primes" << primes.size();
	// cerr << "\nPrimes generated";
	mpz_class t1;

	int trials;
	trials = 3;

	int count, res, size, startIndex, j;
	 j=0;
	 mpz_class* data = genTestData(j);
	 size =200;
	 startIndex = 100;
		 
	 for (int i = startIndex; i < size; ++i)
	 {
		cout << "\n\n\n TRYING TO FACTORIZE NEW NUMBER";
	 	t1 = data[i];
	  	res =factorize(t1, trials);
	  	count = count+res;
	 }

	cout << "\n Could factorize  " << count << " numbers out of 200 starting on number "<<startIndex <<"whith j = " << 0;
	 
	
	

	return (0);
}


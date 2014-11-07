#include <iostream>
#include <stdio.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <limits>
#include <math.h>
#include <map>


using namespace std;


static const  int smallPrimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};
static const int primeLength = 26;
static const int primeSafety = 10;
static const double e = 2.718281;

bool xOr(bool val1, bool val2){
	if(val1 != val2){
		return true;
	}
	else{
		return false;	
	} 
}

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

mpz_class pollard(mpz_class seed, mpz_class add, mpz_class N, int timeOut){
	// check small primes

	std::clock_t pClock;
    double duration;
    pClock = std::clock();

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
	upperBound = sqrt(upperBound);
	bound = (upperBound+lowerBound)/2;
	bound = lowerBound;
		//cerr <<"\n upperBound" <<upperBound;
		//cerr << "\n using bound "<< bound;
	
	product = 1;
	a = seed;
	b = seed;

	
	do {
		product = 1;
		for(double i =0; i<bound;i++){
			a=(pow(a,2) + add) % N;
			b=(pow(b,2) + add) % N;
			b=(pow(b,2) + add) % N;
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
	duration = (std::clock() - pClock ) / (double) CLOCKS_PER_SEC;
		//cerr << "trying...";
	}while(a!=b && duration<timeOut);


	
    //cout<<"Calculation time: "<< duration <<'\n';	



	p = "-1";
	return p;

}
vector<mpz_class> genPrimeBase( int maxB){
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
		//cerr << "q is "<< q<< "\n";
		s=s+1;
	}
	//Find z = a non-quadratic residue in p
	z=2;
	while(legendreSymbol(z,p)!=-1){
		z=z+1;
	}
	c = modExponential(z,q,p);
	res = modExponential(n, (q+1)/2, p);
	t = modExponential(n,q,p);
	m=s;

	while(true){
		//cerr <<"\nInside infinity";
		//cerr << "\n t: " <<t <<" c: "<< c << " res " << res <<" m: "<< m<<" z: "<< z <<" q: "<< q;
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
		//cerr <<"\n after for loop is i: "<<i <<"is m "<<m <<" is t "<< t;

		two =2;

		unsigned long expo = m.get_ui();
		expo = expo-1-i; 
		//cerr <<"\n expo is " << expo;

		b= modExponential(c, pow(two,expo) , p);
		res= (res*b)  % p;

		t = (t*b*b) % p;
		c= (b*b)%p;
		m = i;
		//cerr << "\n b: " <<b <<" c: "<< c << " res " << res <<" m: "<< m<<" t: "<< t;
		
	}
}

mpz_class getPrimeFactor(mpz_class N, int trials, int timeOut){
	
	gmp_randclass rand (gmp_randinit_default);
	mpz_class seed, add, factor, temp,error;
	error = -1;
	bool primFound = false;
	factor = N;
	if(isPrime(factor,10)){
		return factor;
	}

	while(!primFound){
		//cerr << "\n Entering primFound loop with factor: \n " << factor << "\n Trial nr: ";
		for(int i = 0; i < trials; ++i){
			seed = rand.get_z_range(factor-1) +1;
			add = 1; //TODO Hard CODED
			//cerr << i << " ";
			//cerr << "\n Seed is \n " << seed;
			temp = pollard(seed, add, factor,timeOut);
				if(temp > 0){
					factor = temp;
					//cerr << "\nIs the factor " << factor << " considered prime? " << isPrime(factor,10);
					if(factor > 0 && isPrime(factor,primeSafety)){
						return factor;
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
				//Return -1 to signal failing factorization
				return error;
			}

	}
	cout << "\n Factorization failed!";
	return error;

	
}


std::vector<mpz_class> factorize(mpz_class N, int trials){
	
	vector<mpz_class> factors;	
	gmp_randclass rand (gmp_randinit_default);
	
	mpz_class seed, add, factor, temp;
	bool primFound = false;
	while(N!=1){
		factor = N;
		//cerr << "\nFactor in the outer loop is now " << factor;		
		if(isPrime(factor,10)){
			factors.push_back(factor);
			//cerr << "\n Added factor " << factor;
			break;
		}

		while(!primFound){
			//cerr << "\n Entering primFound loop with factor: \n " << factor << "\n Trial nr: ";
			for(int i = 0; i < trials; ++i){
				seed = rand.get_z_range(factor-1) +1;
				add = 1; //TODO Hard CODED
				//cerr << i << " ";
				//cerr << "\n Seed is \n " << seed;
				temp = pollard(seed, add, factor, 10);
					if(temp > 0){
						factor = temp;
						//cerr << "\nIs the factor " << factor << " considered prime? " << isPrime(factor,10);
						if(factor > 0 && isPrime(factor,primeSafety)){
							factors.push_back(factor);
							//cerr << "\n Added factor " << factor;
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
					//Add -1 to end of factors to signal uncomplete factorization
					mpz_class error = -1;
					factors.push_back(error);
					return factors;
				}

			}
			N = N / factor;
			primFound = false;

		}
	return factors;
	}



mpz_class calculatePolynom(mpz_class N, mpz_class sqrtN, mpz_class x){
	mpz_class qx = pow(sqrtN+x, 2) - N;
	return qx;
}
string DecToBin2(int number)
{
    string result = "";

    do
    {
        if ( (number & 1) == 0 )
            result += "0";
        else
            result += "1";

        number >>= 1;
    } while ( number );

    reverse(result.begin(), result.end());
    return result;
}

mpz_class calcNonTrivalFactor(mpz_class N, map<mpz_class, vector<bool> > &nullspace, vector<mpz_class> &smoothPolynoms, vector<mpz_class> &factorBase){
	map<mpz_class, bool>  pivotVariables;
	vector<mpz_class> freeVariables;
	mpz_class res;
	int nSmoothPolynoms = smoothPolynoms.size();
	//Init all varaibles to be free
	for(unsigned i; i<smoothPolynoms.size();++i){
		pivotVariables[smoothPolynoms.at(i)] = false;

	}
	//mark privot variables as false;
	for(unsigned i; i< factorBase.size();++i){
		for (unsigned j = 0; j < smoothPolynoms.size(); ++j)
		{
			if(nullspace[factorBase.at(i)].at(j)){
				//Pivot one found
				pivotVariables[smoothPolynoms.at(j)] = true;
				break;
			}
		}
	}
	//caluclate number of free variables
	
	for(unsigned i; i< smoothPolynoms.size();++i){
		if(!pivotVariables[smoothPolynoms.at(i)]){
			freeVariables.push_back(smoothPolynoms.at(i));
		}
	}
	
	int nFreeVariables = freeVariables.size();
	cerr << "\n The number of free variables are " << nFreeVariables;
	int nOfSoultions = pow(2, nFreeVariables);
	int solutionUnderTest = nOfSoultions -1;

	map<mpz_class,bool> testSolution;

	//Gen and test the differen soulutions
	for(int i =0; i< nOfSoultions;i++){
		//Init testsolution to all ones
		for(int j=0;j<nSmoothPolynoms;++j){
			testSolution[smoothPolynoms.at(j)] = true;
		}

		//SetValue of free variables
		string solutionString = DecToBin2(solutionUnderTest);	
		unsigned nRestingZeros = nFreeVariables - solutionString.size();
		unsigned index = 0;
		mpz_class freeVariable;
		for(index= 0; index< nRestingZeros;index++){
				freeVariable = freeVariables.at(index);
				testSolution[freeVariable] =false;
		}
		
		for(unsigned i = 0; i < solutionString.size(); ++i) {
			if(index>=freeVariables.size()){
				break;
			}

			int a = (int) solutionString[i];
    		if(a == 0){
    			freeVariable = freeVariables.at(index);
    			testSolution[freeVariable] =false;
    		}
    		index++;    		
		}


		//Calculate values for pivot variables
		int initedPivots =0;
		while(initedPivots+nFreeVariables<=nSmoothPolynoms){
			vector<bool> row = nullspace[factorBase.at(initedPivots)];
			bool pivotFound;
			bool pivotValue = false;
			mpz_class pivot;
			for (unsigned k = 0; k < row.size(); ++k)
			{
				if(row.at(k)){
					if(!pivotFound){
						pivotFound = true;
						pivot = smoothPolynoms.at(k);
					}
					else{
						pivotValue = !pivotValue;
					}

				}
			}
			if(pivotFound){
				testSolution[pivot] = pivotValue;
			}

		}
		//Try the generated solution
		mpz_class prod1 =1;
		mpz_class rawPolynom; 
		mpz_class prod2 =1;
		mpz_class diff;
		mpz_class solution;

		for(int i=0; i < nSmoothPolynoms;++i){
			if(testSolution[smoothPolynoms.at(i)]){
				prod1 = prod1*smoothPolynoms.at(i);
				//calc value for sqaure 2
				rawPolynom = smoothPolynoms.at(i)+N;
				rawPolynom = sqrt(rawPolynom);
				prod2 = prod2*rawPolynom;
			}
		}
		mpz_class sqaure;
		sqaure= (prod1-prod2)*(prod1+prod2);
		sqaure = sqaure%N;
		res = gcd(sqaure,N);
		if(res!=1 &&  res!=N){
			return res;
		}
		solutionUnderTest = solutionUnderTest-1;
	}
	res = -1;
	return res;
}

map<mpz_class, vector<bool> > getGausSolutionOfTransposedPoynomMatrix( map<mpz_class, vector<bool> > matrix, vector<mpz_class> &factorBase){
	
	int nPrimesInBase = factorBase.size();
	int nSmoothPolynoms = matrix.size();
	vector<mpz_class> smoothPolynoms;

	for(std::map<mpz_class, vector<bool> >::const_iterator it = matrix.begin(); it != matrix.end(); it++)
	{
		mpz_class key = it->first;
		smoothPolynoms.push_back(key);
	}


	//transpose matrix to transpose
	map<mpz_class, vector<bool> > transpose;
	//prepare one row for each prime in factorBase
	for (int i = 0; i < nPrimesInBase; ++i)
	{
		vector<bool> newRow(nSmoothPolynoms); 
		transpose[factorBase.at(i)] = newRow;

	}
	cerr <<"\n Init transpose matrix with " << transpose.size() <<"rows";
	//Transfer values
	for (int i = 0; i <nPrimesInBase; ++i)
	{
		for (int j = 0; j < nSmoothPolynoms; ++j)
		{
			transpose[factorBase.at(i)].at(j) = matrix[smoothPolynoms.at(j)].at(i);
		}
	}
	cerr << "\n Done Transfered values to transpose.. ";
	cerr << "\n print transpose";
	cerr << "\n The smooth polynoms are";
	for (unsigned i = 0; i < smoothPolynoms.size(); ++i)
	{
		cerr << " " << smoothPolynoms.at(i) << " ";
	}
	cerr<<"\n In transpose is number of rows " <<transpose.size();
	
	for (unsigned i = 0; i < transpose.size(); ++i)
	{	
		cerr << "\n prime p "<< factorBase.at(i) << "have bool vector";
		vector<bool> v = transpose[factorBase.at(i)];
		for (unsigned j = 0; j < v.size(); ++j)
		{
			cerr << " " << v.at(j) << " ";
		}
		cerr <<"\n";

	}


	//Do gaus on transpose;
	
	mpz_class prime; 
	mpz_class polynom;
	//For each column 
	for (int j = 0; j < nSmoothPolynoms; ++j){
		int indexFirstOccurance =-1;
		int indexFirstLeadingOne =-1;
		cerr << "\n try reduce comlumn for poynom " << smoothPolynoms.at(j);
		//for each row
		for (int i = 0; i < nPrimesInBase; ++i){

			//Try find first leading 1 for smoothPolynom j
			bool candidateLeadingOne = true; //sign if this row is a valid candidate for leading 1
			for(int k=0; k<=j;++k){
				 //check that all numbers before transpose[i].at(j) ==0 else is this row no candidate for leading 1
				 if( k<j && (transpose[factorBase.at(i)].at(k)==true) ){
				 		//note the index for first occurance of 1 for smmoth polynom
				 		if(indexFirstOccurance<0){
				 			if(transpose[factorBase.at(i)].at(j)==true){
				 				indexFirstOccurance =i;
				 			}
				 		}
				 	candidateLeadingOne = false;
				 }

				 if(k==j && candidateLeadingOne && transpose[factorBase.at(i)].at(j)){
				 	cerr << "\n First leading 1 for polynom " << smoothPolynoms.at(j) << " is on row with prime " << factorBase.at(i);
				 	indexFirstLeadingOne = i;
				 }
			}

			if(indexFirstLeadingOne>=0){
				//we have fond a leading one
				//break and do gauss;
				break;
			}
		}
		cerr << "\n Time to gauss";
		//Do gauss if leading one found
		if(indexFirstLeadingOne>=0){
			//reduce rows under row with leading on
			for (int i = indexFirstLeadingOne+1; i < nPrimesInBase; ++i)
			{
				prime=factorBase.at(i);
				cerr<< "\n row tested for reduced have prime number " << prime;
				//if row have 1 on possition j do row reduce with help of row that have first leading one 
				if(transpose[factorBase.at(i)].at(j)){
					for(int k=0;k<nSmoothPolynoms;k++){
					polynom = smoothPolynoms.at(k);
						transpose[prime].at(k) = xOr(transpose[factorBase.at(indexFirstLeadingOne)].at(k),transpose[prime].at(k));
						
					}
					cerr<< "\n row with prime number " << prime << " was reduced";

				}
			}
			if(indexFirstOccurance>=0){
			//reduce eventual rows over leading one
			for (int i = indexFirstOccurance; i < indexFirstLeadingOne; ++i)
			{
				//if row have 1 on possition j do row reduce with help of row that have first leading one 
				if(transpose[factorBase.at(i)].at(j)){
					for(int k=0;k<nSmoothPolynoms;k++){
						transpose[factorBase.at(i)].at(k) = xOr(transpose[factorBase.at(indexFirstLeadingOne)].at(k),transpose[factorBase.at(i)].at(k));
					}
				}
			}
			}

		}
		//Continue with the next polynoms

	}
	cerr << "\n The smooth polynoms are";
	for (unsigned i = 0; i < smoothPolynoms.size(); ++i)
	{
		cerr << " " << smoothPolynoms.at(i) << " ";
	}

	cerr<<"\n In transpose is number of rows " <<transpose.size();
	
	for (unsigned i = 0; i < transpose.size(); ++i)
	{	
		cerr << "\n prime p "<< factorBase.at(i) << "have bool vector";
		vector<bool> v = transpose[factorBase.at(i)];
		for (unsigned j = 0; j < v.size(); ++j)
		{
			cerr << " " << v.at(j) << " ";
		}
		cerr <<"\n";

	}

	return transpose;
}





map<mpz_class, vector<bool> > genSmothPolynoms(mpz_class N, mpz_class sqrtN, vector<mpz_class> & factorBase, vector<mpz_class> & x1List, vector<mpz_class> & x2List, vector<double> & logPList){


	std::clock_t pClock;
    double duration;
    pClock = std::clock();


	map<mpz_class, vector<bool> > factoredPolynoms;
	mpz_class xInPolynom =1;
	mpz_class firstXInBatch =1;
	mpz_class p,temp;

	while(factoredPolynoms.size() < factorBase.size()+3 ){ // < (factorBase.size + 10) TODO
	//Create polynoms
	cerr << "\n\n Trying to generate more smooth polynoms \n";
	//int nPolynoms = 21474836;
	int nPolynoms = 100; 
	vector<double> polynomsLog(nPolynoms);

	
	if(firstXInBatch==1){
		for(int i =1; i<=nPolynoms;++i){		
			mpz_class qx = calculatePolynom(N, sqrtN, xInPolynom);
			polynomsLog.at(i-1) = (log(qx.get_d()));
			xInPolynom = xInPolynom+1;
		}
	}
	else{
		mpz_class lastX = firstXInBatch+nPolynoms;
		mpz_class lastQX = calculatePolynom(N,sqrtN, lastX); 
		double lastLogQX = log(lastQX.get_d()); 
		for(int i =1; i<=nPolynoms;++i){		
			polynomsLog.at(i-1) = lastLogQX;
			xInPolynom = xInPolynom+1;
		}

	}
	cerr <<"\n polynoms in batch " << polynomsLog.size();
	//cerr << "\nlowest log polynom " << polynomsLog.front();
	//cerr << "\n highest log polynom " << polynomsLog.back();
	//cerr << "\n diff between lowest and highest polynom " << polynomsLog.back() - polynomsLog.front();
	// - SIEVE TIME -
	//Sieve the polynomsLog
	//Obs ignoring p==2
		
	for(unsigned i = 1; i< factorBase.size(); ++i){
		
		p = factorBase.at(i);

		//Add p until reach x-values in this batch 
		while(x1List.at(i) < firstXInBatch){
			x1List.at(i) = x1List.at(i) +p;
		}
		while(x2List.at(i) < firstXInBatch){
			x2List.at(i) = x2List.at(i) +p;
		}
		
		//cerr << "\n Entering log sieve with p "<< p; 
		mpz_class x1Index = x1List.at(i) - firstXInBatch;
		mpz_class x2Index = x2List.at(i) - firstXInBatch;

		for(unsigned long j = x1Index.get_ui(); j < polynomsLog.size(); j=j+p.get_ui()){

			polynomsLog.at(j) = polynomsLog.at(j) - logPList.at(i);		
		}
		for(unsigned long j = x2Index.get_ui(); j < polynomsLog.size(); j=j+p.get_ui()){
			polynomsLog.at(j) = polynomsLog.at(j) - logPList.at(i);
		}

	}

	double factorLimit = log(factorBase.back().get_d())+1;
	vector<mpz_class> factors;
	
	cerr << "\n done with log sieving factor limit is " << factorLimit;


	// sieve out the smooth polynoms.
	
	int count=0;
	for(unsigned i =0; i<polynomsLog.size();++i){
		//cerr << "\n polynomLog resten är " << polynomsLog.at(i) << " and index is " << i<<" polynom is "<< polynoms.at(i);
		if(polynomsLog.at(i)< factorLimit){
			count++;
			mpz_class originPolynom = calculatePolynom(N,sqrtN,(firstXInBatch+i));
			mpz_class numberToFactor = originPolynom;
			vector<bool> polynomInFactorBase(factorBase.size(), false);	
			while(numberToFactor!= 1){
				temp = getPrimeFactor(numberToFactor,10,20);
				bool validFactor = false;
				//check if recived factor is in factorBase				
				if (temp>0){
					for(unsigned k = 0; k < factorBase.size(); ++k){
						if(factorBase.at(k) == temp){ // Factor in factorBase, prepare vector 
							polynomInFactorBase.at(k) = !polynomInFactorBase.at(k); // flip value
							numberToFactor= numberToFactor/temp;
							validFactor = true;
							break;
						}
					}
				}
				if(validFactor == false){
					break;
				}
			}
			if(numberToFactor ==1){
				factoredPolynoms[originPolynom] = polynomInFactorBase;
			}
		}

	}
		// Factorization of polynoms complete. Time to gauss
		cerr << "\n Number of logsmooth polynoms " << count ;		
		cerr << "\n FactorBase size is " << factorBase.size();
		cerr << "\n Number of polynoms generated " << polynomsLog.size();
		cerr << "\n Number of factored polynoms " << factoredPolynoms.size();
		duration = (std::clock() - pClock ) / (double) CLOCKS_PER_SEC;
		cerr << "\n we have now spent " << duration << " seconds on generating smooth polynoms ";

		firstXInBatch = xInPolynom;
	}
	return factoredPolynoms;


}


int quadraticSieve(mpz_class N, vector<mpz_class> & primes){
	
	mpz_class temp;

	cout << "\n QS trying to factorize \n" << N <<"\n";
	vector<mpz_class> factorBase;
	
	
	
	//Create factorBase
	double B, doubleN, bExp;
	doubleN = N.get_d();
	bExp = sqrt(log(doubleN)*log(log(doubleN)))*0.5;
	B = 3*pow(e, bExp);
	B=43; //TODO remove
	temp = 2;
	factorBase.push_back(temp);
	//cerr <<"\n our factorBase contains ";
	for(unsigned i = 1; i < primes.size(); ++i)
	{	
		if(primes[i]>B){
			break;
		}
		unsigned long exponent = (primes.at(i).get_ui() -1)/2;
		temp = modExponential(N, exponent, primes.at(i));
		if(temp == 1){
			factorBase.push_back(primes.at(i));
			//cerr <<"\n" << primes.at(i);
		} 
	}
	//	Factorbase done
	cerr << "\n Factorbase done, size is " << factorBase.size() << " and largest prime is " << factorBase.back();	
	// precalculate x1,x2 and log(p) for all p's in factorBase
	std::vector<double> logPList(factorBase.size());
	std::vector<mpz_class> x1List(factorBase.size());
	std::vector<mpz_class> x2List(factorBase.size());

	
	// //TODO check if this work
	// mpf_class floatN = N;
	// floatN = sqrt(floatN);
	// floatN = floor(floatN);
	// mpz_class sqrtN(floatN);
	mpz_class sqrtN = sqrt(N);
	sqrtN = floor(sqrtN);
	

	mpz_class R1, R2,p,x1,x2;
	for(unsigned i = 1; i< factorBase.size(); ++i){
		
		p = factorBase.at(i);
		//cerr << "\n Entering log sieve with p "<< p; 
		
		R1 = shanksTonelli(N, p);
		//cerr << "\n done shanksTonelli for p " <<p;
		R2 = p - R1;
		x1 = (R1-sqrtN) % p;
		x2 = (R2 -sqrtN) % p;
		double pLog = log(p.get_d());		

		logPList.at(i) = pLog;
		x1List.at(i) = x1;
		x2List.at(i) = x2;
	}

	map<mpz_class, vector<bool> > smoothPolynoms =  genSmothPolynoms(N,sqrtN,factorBase, x1List, x2List,logPList);
	cerr <<"\n number of generated smoothPolynoms is " << smoothPolynoms.size();
	cerr <<"\n number of primes in factorBase is " << factorBase.size();

	map<mpz_class, vector<bool> > nullspace = getGausSolutionOfTransposedPoynomMatrix(smoothPolynoms,factorBase);
	return 0;
}


void testGetPrimeFactor(){
	mpz_class res;
	mpz_class N = 57385973589534549;
	res = getPrimeFactor(N,10,20);
	cerr <<"\n found factor "<< res << " in " << N;
}
void testGenNullSpace(){
	int maxB = 50000;
	cerr << "\nGenerating primes...";
	vector<mpz_class> primes;	
	primes = genPrimeBase(maxB);
	cerr <<"\nsize of primes " << primes.size();
	cerr << "\nPrimes generated";
	
	mpz_class n = 90283;
	quadraticSieve(n, primes);
	

}


int main() {	
	
	// // //Gen test data
	// int j =0;		
	// mpz_class* data = genTestData(j);

	// //Gen prime base
	// int maxB = 5000000;
	
	// cerr << "\nGenerating primes...";
	// vector<mpz_class> primes;	
	// primes = genPrimeBase(maxB);
	// cerr <<"\nsize of primes " << primes.size();
	// cerr << "\nPrimes generated";
	
	
	// int trials, count, stopIndex, startIndex;
	// trials = 3;
	// stopIndex =12;
	// startIndex = 9;
	// vector<mpz_class> factors;

	//  for (int i = startIndex; i < stopIndex; ++i)
	//  {
	// 	cout << "\n\n\n TRYING TO FACTORIZE NEW NUMBER";
	//  	cout << "\n Trying to factorize \n" <<data[i] <<"\n";
	// 	std::clock_t start;
 //    	double duration;
 //    	start = std::clock();


	//   	factors = factorize(data[i], trials);
	//   	if(factors.back() == -1){
	//   		cerr << "\n Uncomplete factorization ";
	//   	}
	//   	else{
	//   		count++;
	//   	}
	//   	cout << "\n The number of prime factors is " << factors.size(); 
	// 	//print factors
	// 	for (unsigned i = 0; i < factors.size() ; ++i)
	// 	{
	// 		cout <<"\n " << factors.at(i);
	// 	}
	// 	cout << "\n";
	// 	duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//     cout<<"\n Calculation time: "<< duration <<'\n';	
	//  }

	// cout << "\n Could factorize  " << count << " numbers out of " << stopIndex-startIndex << " starting on number "<<startIndex <<" with j = " << j;
	 
	//TEST CODE
	// mpz_class t1;
	// //t1 = "66293490818913051990436158760858275128292159574918701991";
	// t1 = 90283;
	
	// quadraticSieve(t1,primes);
	//testGetPrimeFactor();
	testGenNullSpace();

	return (0);
}



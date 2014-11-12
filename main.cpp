#include <iostream>
#include <stdio.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <limits>
#include <math.h>
#include <map>
#include <bitset>


using namespace std;


static const  int smallPrimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};
static const int primeLength = 26;
static const int primeSafety = 10;
static const double e = 2.718281;



vector<unsigned long> rowXOR(vector<unsigned long> & v1, vector<unsigned long> & v2){
	vector<unsigned long> result;
	for(int i = 0; i < v1.size(); ++i){
		result.push_back(v1.at(i) xor v2.at(i));
	}
	return result;
}

bool getBitAt(vector<unsigned long> & words, unsigned long bitIndex){
	//cerr << "\n Inside getBitAt with length of words" << words.size() << " and bitIndex " << bitIndex; 
	unsigned wordSize = 64;
	unsigned wordIndex = 0;
	while(bitIndex >= 64){
		bitIndex -= 64;
		wordIndex++;
	}
	unsigned long word = words.at(wordIndex);
	return ((word >> bitIndex) & 1);
}

vector<bool> getBoolVector(vector<unsigned long> words, unsigned nBits){
	vector<bool> result;
	for(unsigned i = 0; i < nBits; ++i){
		result.push_back(getBitAt(words, i));

	}
	return result;
}

bitset< 64 > getBitset(vector<bool> & bools, unsigned startIndex, unsigned stopIndex){
	bitset< 64 > result;
	for(unsigned i = startIndex; i < stopIndex; ++i ){
		if(bools[i]){
			result[i] = 1;
		} else {
			result[i] = 0;
		}
	}
	return result;
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

mpz_class getPrimeFactorByPollard(mpz_class N, int trials, int timeOut){
	
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
			cout << "\n Factorization failed by Pollard!";
			//Return -1 to signal failing factorization
			return error;
		}

	}
	cout << "\n Factorization failed by Pollard!";
	return error;

	
}

mpz_class calculatePolynom(mpz_class N, mpz_class sqrtN, mpz_class x){
	mpz_class qx = pow(sqrtN+x, 2) - N;
	return qx;
}
string DecToBin2(int number){
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
	unsigned nSmoothPolynoms = smoothPolynoms.size();
	vector<mpz_class> rowsWithPivot;

	//cerr << "\n entering calc non trival factor with a factorBase with size " << factorBase.size() << " smoothPolynoms size " << smoothPolynoms.size();
	//Init all varaibles to be free
	for(unsigned i; i<smoothPolynoms.size();++i){
		pivotVariables[smoothPolynoms.at(i)] = false;
	}
	//cerr << "\n size of variables are " << pivotVariables.size() ;
	//mark privot variables as not free;
	for(unsigned i=0; i< factorBase.size();++i){
		for (unsigned j = 0; j < smoothPolynoms.size(); ++j)
		{
			if(nullspace[factorBase.at(i)].at(j)){
				//Pivot one found
				//cerr << "\n" << smoothPolynoms.at(j) << " was found to be pivot";
				rowsWithPivot.push_back(smoothPolynoms.at(j));
				pivotVariables[smoothPolynoms.at(j)] = true;
				break;
			}
		}
	}
	//caluclate number of free variables
	
	for(unsigned i=0; i< smoothPolynoms.size();++i){
		if(!pivotVariables[smoothPolynoms.at(i)]){
			freeVariables.push_back(smoothPolynoms.at(i));
		}
	}
	
	unsigned nFreeVariables = freeVariables.size();
	cerr << "\n The number of free variables are " << nFreeVariables;
	mpz_class nOfSoultions = pow(2, nFreeVariables);
	mpz_class solutionUnderTest = nOfSoultions -1;
	cerr << "\n number of solutions to test is " << solutionUnderTest;
	map<mpz_class,bool> testSolution;
	
	//Gen and test the differen soulutions
	while(solutionUnderTest>0){
		//cerr << "\n solutionUnderTest is " << solutionUnderTest;
		//Init testsolution to all ones
		for(unsigned j=0;j<nSmoothPolynoms;++j){
			testSolution[smoothPolynoms.at(j)] = true;
		}

		//SetValue of free variables
		//string solutionString = DecToBin2(solutionUnderTest);
		string solutionString =  solutionUnderTest.get_str(2);	
		cerr << "\n bin representation is  " << solutionString;
		
		unsigned nRestingZeros = nFreeVariables - solutionString.size();
		unsigned index = 0;
		mpz_class freeVariable;
		for(index= 0; index< nRestingZeros;index++){
				freeVariable = freeVariables.at(index);
				testSolution[freeVariable] =false;
		}
		
		for(unsigned i = 0; i < solutionString.size(); ++i) {
			if(index>=freeVariables.size()){
				//cerr << " index greater then number of free variables";
				break;
			}

			char a = solutionString[i];
			//cerr << " a is " << a;
    		if(a == '0'){
    			freeVariable = freeVariables.at(index);
    			testSolution[freeVariable] =false;
    		}
    		index++;    		
		}


		//cerr << "\n calc value for pivot variables ";
		//Calculate values for pivot variables
		for(unsigned m =0; m<rowsWithPivot.size();++m){
			vector<bool> row = nullspace[rowsWithPivot.at(m)];
			mpz_class pivot;
			bool pivotFound =false;
			bool pivotValue;
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
		//cerr << "\n generated solution gave the polynoms following value";

		// for (unsigned k = 0; k < smoothPolynoms.size(); ++k)
		// {
		// 	cerr << "\n The Polynom " << smoothPolynoms.at(k) << " got value " << testSolution[smoothPolynoms.at(k)];
		// }

		mpz_class prod1 =1;
		mpz_class rawPolynom; 
		mpz_class prod2 =1;
		mpz_class diff;
		mpz_class solution;

		for(unsigned i=0; i < nSmoothPolynoms;++i){
			if(testSolution[smoothPolynoms.at(i)]){
				prod1 = prod1*smoothPolynoms.at(i);
				//calc value for sqaure 2
				rawPolynom = smoothPolynoms.at(i)+N;
				rawPolynom = sqrt(rawPolynom);
				prod2 = prod2*rawPolynom;
			}
		}
		//cerr << "\n Prod 1 is " << prod1 << " prod2 is " << prod2;
		mpz_class sqaure;
		sqaure= (prod1-prod2)*(prod1+prod2);
		sqaure = sqaure;
		res = gcd(sqaure,N);
		cerr << "\n CALCULATED RES IS " << res;

		if(res!=1 &&  res!=N){
			return res;
		}
		solutionUnderTest = solutionUnderTest-1;
	}
	res = -1;
	return res;
}

mpz_class matrixToFactor( mpz_class N, map<mpz_class, vector<bool> > &matrix, vector<mpz_class> &factorBase){
	
	int nPrimesInBase = factorBase.size();
	int nSmoothPolynoms = matrix.size();
	vector<mpz_class> smoothPolynoms;

	for(std::map<mpz_class, vector<bool> >::const_iterator it = matrix.begin(); it != matrix.end(); it++)
	{
		mpz_class key = it->first;
		smoothPolynoms.push_back(key);
	}

	cerr << "\n Starting to transpose matrix ...";
	//transpose matrix to transpose
	map<mpz_class, vector<bool> > transpose;
	//prepare one row for each prime in factorBase
	for (int i = 0; i < nPrimesInBase; ++i)
	{
		vector<bool> newRow(nSmoothPolynoms); 
		transpose[factorBase.at(i)] = newRow;

	}
	//cerr <<"\n Init transpose matrix with " << transpose.size() <<"rows";
	//Transfer values
	for (int i = 0; i <nPrimesInBase; ++i)
	{
		for (int j = 0; j < nSmoothPolynoms; ++j)
		{
			transpose[factorBase.at(i)].at(j) = matrix[smoothPolynoms.at(j)].at(i);
		}
	}
	cerr << "\n Done Transfered values to transpose ... ";
	map<mpz_class, vector<unsigned long> > transposeWithWords;
	// Take transpose and translate each row from bool vectors to unsigned longs
	// Each prime represents one row
	cerr << "\nTranslating transpose into word 64 based matrix, with nSmoothPolynoms = " << nSmoothPolynoms;
	for(int i = 0; i < nPrimesInBase; ++i){
		vector<unsigned long> wordRow;
		for (int j = 0; j < nSmoothPolynoms; j=j+64){ // enough words to cover all + a not finished word
			bitset< 64 > bits = getBitset(transpose[i], j, j+64);
			wordRow.push_back(bits.to_ulong());
		}
		transposeWithWords[i] = wordRow;
	}
	cerr << "\n Done with translation ...";








	//Do gaus on transpose;
	cerr << "\n Start gauss";
	mpz_class prime; 
	mpz_class polynom;
	//For each column 
	cerr << "\n number of rows (primes) are " << factorBase.size();
	cerr << "\n number of rows in transposeWithWords " << transposeWithWords.size();
	cerr << "\n number of columns (polynoms) are " << nSmoothPolynoms;
	cerr << "\n Now done with column :";

	for (int j = 0; j < nSmoothPolynoms; ++j){ // for each column

		if(j % 100 == 0){
		cerr << " " << j;
		}

		int indexFirstOccurance =-1;
		int indexFirstLeadingOne =-1;
		//cerr << "\n try reduce column for polynom " << smoothPolynoms.at(j);
		//for each row
		for (int i = 0; i < nPrimesInBase; ++i){

			//Try finding first leading 1 for column smoothPolynom j
			bool candidateLeadingOne = true; //sign if this row is a valid candidate for leading 1
			for(int k=0; k<=j;++k){
				 //check that all numbers before transpose[i].at(j) ==0 else this row is no candidate for leading 1
				 if( k<j && (getBitAt(transposeWithWords[i], k)==true) ){
				 		//note the index for first occurance of 1 for smooth polynom
				 		if(indexFirstOccurance<0){
				 			if(getBitAt(transposeWithWords[i],j)==true){
				 				indexFirstOccurance =i;
				 			}
				 		}
				 	candidateLeadingOne = false;
				 }
				 if(k==j && candidateLeadingOne && getBitAt(transposeWithWords[i],j)){
				 	//cerr << "\n First leading 1 for polynom " << smoothPolynoms.at(j) << " is on row with prime " << factorBase.at(i);
				 	indexFirstLeadingOne = i;
				 }
			}

			if(indexFirstLeadingOne>=0){
				//we have found a leading one
				//break and do gauss;
				break;
			}
		}
		//cerr << "\n Time to gauss";
		//Do gauss if leading one found
		if(indexFirstLeadingOne>=0){
			//reduce rows under row with leading on
			for (int i = indexFirstLeadingOne+1; i < nPrimesInBase; ++i) // for each row
			{
				//cerr<< "\n row tested for reduced have prime number " << prime;
				//if row have 1 on position j do row reduce with help of row that have first leading one 
				if(getBitAt(transposeWithWords[i],j)){
					//perform XOR between the row with leading 1 and the one which also has 1
					transposeWithWords[i] = rowXOR(transposeWithWords[indexFirstLeadingOne], transposeWithWords[i]);
					//cerr<< "\n row with prime number " << prime << " was reduced";

				}
			}
			if(indexFirstOccurance>=0){
			//reduce eventual rows over leading one
			for (int i = indexFirstOccurance; i < indexFirstLeadingOne; ++i){
				//if row have 1 on possition j do row reduce with help of row that have first leading one 
				if(getBitAt(transposeWithWords[i],j)){
					transposeWithWords[i] = rowXOR(transposeWithWords[indexFirstLeadingOne], transposeWithWords[i]);
					}
				}
			}

		}
		//Continue with the next column (polynom)

	}
	cerr << "\n Done gaussing, translating back to bool vector from word vector ... ";
	// for each row in transposeWithWords translate its values to transpose (overwriting old values pre gauss)
	for(int i = 0; i < factorBase.size(); ++i){
		transpose[i] = getBoolVector(transposeWithWords[i], nSmoothPolynoms);
	}
	cerr << "\n Done translating back, ..."; 


	mpz_class res = calcNonTrivalFactor(N, transpose,smoothPolynoms,factorBase);

	return res;
}





map<mpz_class, vector<bool> > genSmothPolynoms(mpz_class N, mpz_class sqrtN, vector<mpz_class> & factorBase, vector<mpz_class> & x1List, vector<mpz_class> & x2List, vector<double> & logPList){


	std::clock_t pClock;
    double duration;
    pClock = std::clock();


	map<mpz_class, vector<bool> > factoredPolynoms;
	mpz_class xInPolynom =1;
	mpz_class firstXInBatch =1;
	mpz_class p,temp;

	while(factoredPolynoms.size() < factorBase.size()+10){
	//Create polynoms
	//cerr << "\n\n Trying to generate more smooth polynoms \n";
	int nPolynoms = factorBase.size()*200;//: 21474836;
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
	//cerr <<"\n polynoms in batch " << polynomsLog.size();
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

	double factorLimit = log(factorBase.back().get_d())+1; //TODO
	
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
				temp = getPrimeFactorByPollard(numberToFactor,10,20);
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


mpz_class quadraticSieve(mpz_class N, vector<mpz_class> &primes){
	
	mpz_class temp;

	cout << "\n QS trying to factorize \n" << N <<"\n";
	vector<mpz_class> factorBase;
	
	
	
	//Create factorBase
	double B, doubleN, bExp;
	doubleN = N.get_d();
	bExp = sqrt(log(doubleN)*log(log(doubleN)))*0.5;
	B = 2*pow(e, bExp);
	//TODO
	B=43;
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
	//qrtN = floor(sqrtN);
	

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

	cerr << "\n \n Smooth polynom generation done, starting with gaussian elimination for number \n " << N;
	
	mpz_class res = matrixToFactor(N, smoothPolynoms , factorBase);
	return res;
}


std::vector<mpz_class> factorize(mpz_class N, vector<mpz_class> &primes ){
	cout << "\n Trying to factorize \n" <<N <<"\n";		
	vector<mpz_class> factors;	
	mpz_class temp,factor;
	int trials = 2;
	int timeOut = 60;

	//Check if N is prime
	if(isPrime(N,10)){
		factors.push_back(N);
		return factors;
	}
	//Do pollard to find small factors
	temp =1;
	while(temp>0){
		temp = getPrimeFactorByPollard(N, trials, timeOut);
		if(temp>0){
			factor=temp;
			cerr <<"\n Found factor " << temp << " by pollard";
			N = N / factor;
			factors.push_back(factor);
			cerr << "\n N is now " << N;
		}
		if(N==1){
			return factors;
		}
	}
	//check if result is 1 or peime
	if(isPrime(N,10)){
		factors.push_back(N);
		return factors;
	}
	//try quadratic sive on rest;
	
	while(N!=1){
		factor = quadraticSieve(N,primes);
		if(factor<0){
			cerr <<"\n Found factor " << temp << " by quadratic sieve";
			N = N / factor;
			factors.push_back(factor);
			cerr << "\n N is now" << N;
		}
		else{
			cerr <<"\n  Quadratic sieve failed to factor " << N;
			factor = -1;
			factors.push_back(factor);
		}
		if(isPrime(N,10)){
			factors.push_back(N);
			return factors;
		}
	}
	return factors;
	}

	void testQuadraticSieveOn90283() {
	int maxB = 50000;
	cerr << "\nGenerating primes...";
	vector<mpz_class> primes;	
	primes = genPrimeBase(maxB);
	cerr <<"\nsize of primes " << primes.size();
	cerr << "\nPrimes generated";
	
	mpz_class n = 90283;
	mpz_class res;
	res = quadraticSieve(n, primes);
	cerr <<"\n QS found factor " << res << " in " << n;

	}


// void testGetPrimeFactor(){
// 	mpz_class res;
// 	mpz_class N = "57385973589534549";
// 	res = getPrimeFactorByPollard(N,10,20);
// 	cerr <<"\n found factor "<< res << " in " << N;
// }
// void testGenNullSpace(){
// 	int maxB = 50000;
// 	cerr << "\nGenerating primes...";
// 	vector<mpz_class> primes;	
// 	primes = genPrimeBase(maxB);
// 	cerr <<"\nsize of primes " << primes.size();
// 	cerr << "\nPrimes generated";
	
// 	mpz_class n = 90283;
// 	quadraticSieve(n, primes);
	

// }

void testGetBinaryConversion(){
	std::vector<bool> b;
	b.push_back(true);
	b.push_back(false);
	b.push_back(false);
	b.push_back(true);
	b.push_back(true);
	unsigned startIndex = 0;
	unsigned stopIndex = b.size();
	bitset< 64 > output;
	output = getBitset(b, startIndex, stopIndex);
	cout << "\n bitset printed is " << output.to_string();
	unsigned long res = output.to_ulong();
	cout << "\n unsigned long is printed " << res;
	cout << "\n bit number 0 is " << ((res >> 0) & 1); 
	cout << "\n bit number 1 is " << ((res >> 1) & 1); 
	cout << "\n bit number 2 is " << ((res >> 2) & 1); 
	cout << "\n bit number 3 is " << ((res >> 3) & 1); 
	cout << "\n bit number 4 is " << ((res >> 4) & 1); 


}

void testGetBitAt(){
	bitset< 64 > b1(14);
	bitset< 64 > b2(1025);
	cerr << "\n First word is " << b1.to_string();
	cout << "\n Second word is " << b2.to_string();
	vector<unsigned long> v;
	v.push_back(b1.to_ulong());
	v.push_back(b2.to_ulong());
	cerr << "\n bit with index 0 is " << getBitAt(v, 0);
	cerr << "\n bit with index 1 is " << getBitAt(v, 1);
	cerr << "\n bit with index 64 is " << getBitAt(v, 64);
	cerr << "\n bit with index 65 is " << getBitAt(v, 65);

}

void testUnsignedLongSize(){
	cerr << "size of unsigned long is " << sizeof(unsigned long)*8;
}

int main() {	

	testQuadraticSieveOn90283();
	
	testGetBitAt();

	//Gen test data
		int j =0;		
	mpz_class* data = genTestData(j);

	// Test binary conversion methods:
	//testGetBinaryConversion();
	//testGetBitAt();



	//Gen prime base
	int maxB = 5000000; //TODO add 2 0
	cerr << "\nGenerating primes...";
	vector<mpz_class> primes;	
	primes = genPrimeBase(maxB);
	cerr <<"\nsize of primes " << primes.size();
	cerr << "\nPrimes generated";

	
	int count, stopIndex, startIndex;
	//trials = 3;
	startIndex = 25;
	stopIndex =28;
	
	vector<mpz_class> factors;

	 for (int i = startIndex; i < stopIndex; ++i)
	 {
		cout << "\n\n\n TRYING TO FACTORIZE NEW NUMBER";
	 	std::clock_t start;
    	double duration;
    	start = std::clock();


	  	factors = factorize(data[i], primes);
	  	if(factors.back() == -1){
	  		cerr << "\n Uncomplete factorization ";
	  	}
	  	else{
	  		count++;
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
	 }

	cout << "\n Could factorize  " << count << " numbers out of " << stopIndex-startIndex << " starting on number "<<startIndex <<" with j = " << j;
	 
	// //TEST CODE

	//mpz_class t1;
	//t1 = "90283";	
	//quadraticSieve(t1,primes);
	//testGetPrimeFactor();
	// testGenNullSpace();

	return (0);
}



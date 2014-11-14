#include <iostream>
#include <stdio.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <limits>
#include <math.h>
#include <map>
#include <bitset>
#include <list>


using namespace std;


static const  int smallPrimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};
static const int primeLength = 26;
static const int primeSafety = 10;
static const double e = 2.718281;
static const int bitsOfUL = sizeof(unsigned long) *8;
static const int bytesOfUL = sizeof(unsigned long);

bitset<bitsOfUL > getBitset(vector<bool> & bools, unsigned startIndex, unsigned stopIndex){
	bitset< bitsOfUL > result;
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


class BitMatrix
{
	
	public:
		BitMatrix(mpz_class numberToFactor, map<mpz_class,std::vector<bool> > &rawMatrix, int sizeOfFactorBase);
		
		~BitMatrix();
		void printPrimeMatrix();
		void printPolyMatrix();
		void printZeros();
		mpz_class getFactorByQS();
		void gaussBitMatrix();
		void setZeros();

	private:
		mpz_class N;
		int nPrimesInBase;
		bitset<bitsOfUL> **bitMatrix;
		int nSmoothPolynoms;
		int nPrimeCols;
		int nPolyCols;
		int nCols;
		vector<mpz_class> smoothPolynoms;
		vector<int> zeros;
};

mpz_class BitMatrix::getFactorByQS(){
	printPrimeMatrix();
	gaussBitMatrix();
	cerr << "\n  After gauss";
	printPrimeMatrix();
	printPolyMatrix();
	setZeros();
	mpz_class res=-1;
	mpz_class square;
	vector<mpz_class> polynomsInSolution;
	mpz_class prod1 =1;
	mpz_class rawPolynom; 
	mpz_class prod2 =1;
	mpz_class diff;

	for(unsigned i=0; i<zeros.size();++i){
		int col = nPrimeCols;
		prod1 =1;
		prod2 =1;
		//find polynoms in soulution		
		for(int j=0;j<nSmoothPolynoms;++j){
			if(j%bitsOfUL){
				col++;
			}
			if(bitMatrix[zeros.at(i)][col].test(j%bitsOfUL)){
				polynomsInSolution.push_back(smoothPolynoms.at(j));
			}			
		}
		//try soulution
		for(unsigned j=0;j<polynomsInSolution.size();++j){
			prod1 = prod1*polynomsInSolution.at(j);
			rawPolynom = polynomsInSolution.at(j)+N;
			rawPolynom = sqrt(rawPolynom);
			prod2 = prod2*rawPolynom;
		}
		square = (prod1-prod2)*(prod1+prod2);
		res = gcd(square,N);
		if(res!=1 &&  res!=N){
			cerr << "\n CALCULATED RES IS " << res;
			return res;
		}

	}
	cerr << "\n FACTORIZE FAILED BY QS ";
	res = -1; 
	return res;
}

void BitMatrix::setZeros(){
	
	bool foundZero;
	for(int i =0;i < nSmoothPolynoms; i++){
		foundZero=true;
		int nBitsSet =0;
		for(int j=0;j<nPrimeCols;j++){
			if(bitMatrix[i][j].count()>0){
				foundZero = false;
				break;
			}				
		}
		if(foundZero){
			zeros.push_back(i);
		}
	}
}
//GaussMethod

void BitMatrix::gaussBitMatrix(){
	
	bool possibleLeadingOneRows[nSmoothPolynoms];

	//Create linked list with all row index
	for (int i = 0; i < nSmoothPolynoms; ++i)
		possibleLeadingOneRows[i] = true;	
	
	int indexLeadingOne=0;
	int colToGauss=0;
	int bitToGauss;
	//Find leading one
	for(int i=0; i<nPrimesInBase;i++){
		bitToGauss = i ;
		bool foundLeadingOne=false;
		
		indexLeadingOne =0;

		if(bitToGauss>0&&bitToGauss%bitsOfUL==0)
			colToGauss++;
		//Find leading one
		for(int j=0;j<nSmoothPolynoms;++j){
			if(possibleLeadingOneRows[indexLeadingOne] && bitMatrix[indexLeadingOne][colToGauss].test(bitToGauss%bitsOfUL)){
				possibleLeadingOneRows[indexLeadingOne] =false;
				foundLeadingOne = true;
				break;
			}
			indexLeadingOne++;
		}
		if(foundLeadingOne)
			for(int j=0;j<nSmoothPolynoms;++j)
				if(j!=indexLeadingOne)
					if(bitMatrix[j][colToGauss].test(bitToGauss%bitsOfUL))
						for(int k= colToGauss;k<nCols;k++)
							bitMatrix[j][k] ^= bitMatrix[indexLeadingOne][k];							

	}	 
}

//Print methods
void BitMatrix::printPrimeMatrix(){
	cerr <<"\n Print Prime matrix \n";
	for (int i = 0; i < nSmoothPolynoms; ++i)
	{
		cerr <<"\n Row " << i << ": ";
		for(int j=0;j<nPrimeCols;++j)
			cerr << bitMatrix[i][j].to_string() << " ";
	}

}

void BitMatrix::printPolyMatrix(){
	cerr <<"\n Print Polynom matrix \n";
	for (int i = 0; i < nSmoothPolynoms; ++i)
	{
		cerr <<"\n Row " << i << ": ";
		for(int j=nPrimeCols;j<nCols;++j)
			cerr << bitMatrix[i][j].to_string() << " ";
	}

}
void BitMatrix::printZeros(){
	cerr <<"\n Row with zeros have index";		
	for(unsigned i=0;i<zeros.size();++i){
		cerr <<  zeros.at(i) << " ";
	}
}


//Constructor
BitMatrix::BitMatrix(mpz_class numberToFactor, map<mpz_class,std::vector<bool> > &rawMatrix,int sizeOfFactorBase){
	N=numberToFactor;
	nPrimesInBase = sizeOfFactorBase;
	nSmoothPolynoms = rawMatrix.size();
	nPrimeCols = ceil((double) nPrimesInBase/(double) bitsOfUL);

	nPolyCols = ceil((double) nSmoothPolynoms/(double) bitsOfUL);
	nCols = nPrimeCols+nPolyCols;
	//cerr << "\n Bits of ul is " << bitsOfUL;
	//cerr << "\n nPrimesInBase is " << nPrimesInBase << " nSmoothPolynoms: " << nSmoothPolynoms << " nPrimeCols " << nPrimeCols  <<" nPolyCols "<< nPolyCols <<"\n"; 
	//Declare bitMatrix
	bitMatrix = new bitset<bitsOfUL> *[nSmoothPolynoms];

	for (int i = 0; i < nSmoothPolynoms; ++i)
		bitMatrix[i] = new bitset<bitsOfUL>[nCols];

	//fill with data
	vector<bool> rawRow;
	int k;		
	int i =0;
	
	int polyColPos = nCols;
		
	//cerr << "\n Before loop";
			
	for(std::map<mpz_class, vector<bool> >::const_iterator it = rawMatrix.begin(); it != rawMatrix.end(); it++)
	{
		//cerr << "\n loop " << i;
		//Move bools from rawMatrix
		mpz_class key = it->first;	
		smoothPolynoms.push_back(key);
		rawRow = it->second;

		for (int j = 0; j < nPrimeCols; ++j)
		{
			bitset< bitsOfUL > bits = getBitset(rawRow, (unsigned) j*bitsOfUL, (unsigned) j*bitsOfUL+bitsOfUL);
			//cerr << "\n got bitset " << bits.to_string();
			bitMatrix[i][j] = bits; 
		}
		//cerr << "\n Succesful moved rawMatrix row "<< i;
	
 			
 		if(i%bitsOfUL==0){
 			polyColPos = polyColPos -1;
 		}
 		bitset<bitsOfUL> polyBitSet; 
		polyBitSet.set(i%bitsOfUL, true);		
		bitMatrix[i][polyColPos] = polyBitSet;
		i++;

		//cerr << "\n Succesful added Polyending row "<< i;
	
	}
	cerr << "\n construction of BitMatrix done ";
}
//Destructor

BitMatrix::~BitMatrix(){
	
	for (int i = 0; i < nSmoothPolynoms; ++i)
		delete [] bitMatrix[i];

	delete [] bitMatrix;
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


	mpz_class sqrtN = sqrt(N);
	
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
	
	BitMatrix bm(N,smoothPolynoms, (int) factorBase.size());
	mpz_class res = bm.getFactorByQS();
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


void testUnsignedLongSize(){
	cerr << "size of unsigned long is " << sizeof(unsigned long)*8;
}




void testConstructBitMatrix(){
		map<mpz_class,std::vector<bool> > testRawMatrix;
		mpz_class N = 4711;
		
		for(int i =0;i<65;i++){
			mpz_class key = i;
			std::vector<bool> value(4);
			for(int j=0;j<4;j++){
				if(j%2!=0){
					value[j]=true;
				}
				else{
					value[j] =false;
				}
			}
			testRawMatrix[key] = value;
		}

		BitMatrix bitMatrix(N,testRawMatrix,4);
		bitMatrix.printPrimeMatrix();
		bitMatrix.printPolyMatrix();
		bitMatrix.gaussBitMatrix();
		bitMatrix.printPrimeMatrix();
		bitMatrix.printPolyMatrix();
		bitMatrix.setZeros();
		bitMatrix.printZeros();
}

int main() {	


	//testConstructBitMatrix();
	//exit(0);
	//testSetBit();
	testQuadraticSieveOn90283();
	
	//testGetBitAt();

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



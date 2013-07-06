#include <iostream>
#include <fstream>
#include <math.h>
#include <utility>
#include <time.h>

using namespace std;

const unsigned short MAX_POWER = 30;

//returns the value of b^e
unsigned long long mypow(unsigned int b, unsigned int e)
{
  unsigned long long ans = 1;
  for (unsigned int i=0;i<e;i++)
    ans *= b;
  return ans;
}

//determines the number of digits needed to represent n in binary
unsigned short getNumBitsNeeded(unsigned long long n)
{
  unsigned int ans=1;
  while (mypow(2,ans) <= n)
    ans++;
  return ans;
}


//returns the log base 2 of x
double lg2(int x)
{
  return log(x)/log(2);
}

//returns the maximum i such that 2^i < x
int deg(int x)
{
  int maxPower =0;
  while (pow(2,maxPower) < x)
    maxPower++;
  return maxPower-1;
}

//returns the negligibility of n by table lookup
unsigned long long int getNegFromTable(unsigned long long n)
{ 
  FILE* negTable;
  negTable = fopen("negTable.bin", "rb");
  unsigned long long ans;

  fseek(negTable, n*sizeof(unsigned long long), SEEK_SET);
  fread(&ans,1,sizeof(unsigned long long), negTable);
  fclose(negTable);
  return ans;
}

//returns a pointer to an array.
//there are maxp_c elements in the array.
//the i'th element represents the number of integers n such that p(n) < i
unsigned long long* getLessThanArray(int maxp_c)
{
  unsigned long long* ans;
  ans = new unsigned long long[maxp_c];
  for (int i=0;i<maxp_c;i++)
    ans[i]=0;

  FILE* negTable;
  negTable = fopen("negTable.bin", "rb");
  unsigned int ARR_SIZE = pow(2,15);
  cout << "ARR_SIZE IS " << ARR_SIZE << endl;
  unsigned long long* arr;
  arr = new unsigned long long[ARR_SIZE];
  unsigned int buffSize = ARR_SIZE * sizeof(unsigned long long);


  unsigned int stop = pow(2,maxp_c - 15);
  for (unsigned int i=0; i < stop; i++)
  {
    if (i!=0) 
      fseek(negTable, ARR_SIZE, SEEK_CUR);

    fread(arr,1, buffSize, negTable);
    for (unsigned long long j=0;j<ARR_SIZE;j++)
      for (int k=0;k<maxp_c;k++)
	if (arr[j] < k)
	  ans[k]++;

  }
  delete[] arr;
  fclose(negTable);
  return ans;  
}

//returns a pointer to an array.
//there are maxp_c elements in the array.
//the i'th element represents the number of integers n such that p(n) = i
unsigned long long* getEqualsToArray(int maxp_c)
{
  unsigned long long* ans;
  ans = new unsigned long long[maxp_c];
  for (int i=0;i<maxp_c;i++)
    ans[i]=0;

  FILE* negTable;
  negTable = fopen("negTable.bin", "rb");
  unsigned int ARR_SIZE = pow(2,15);
  cout << "ARR_SIZE IS " << ARR_SIZE << endl;
  unsigned long long* arr;
  arr = new unsigned long long[ARR_SIZE];
  unsigned int buffSize = ARR_SIZE * sizeof(unsigned long long);

  unsigned int stop = pow(2,maxp_c - 15);
  for (unsigned int i=0; i < stop; i++)
  {
    if (i!=0) 
      fseek(negTable, ARR_SIZE, SEEK_CUR);

    fread(arr,1, buffSize, negTable);
    for (unsigned long long j=0;j<ARR_SIZE;j++)
      for (int k=0;k<maxp_c;k++)
	if (arr[j] == k)
	  ans[k]++;

  }
  delete[] arr;
  fclose(negTable);
  return ans;  
}

//returns the location of the binary representation of n in a lookup table
//the location is given by...
// 1 + sum_{i=1}^{numBits-1} i*2^{i-1}
//^That, plus numBits * n-2^{numBits}
unsigned long long whereToLook(unsigned long long n)
{
  unsigned short numBits = getNumBitsNeeded(n);
  unsigned long long ans = mypow(2,numBits-1)*numBits - mypow(2,numBits) + 2 + numBits * (n-mypow(2,numBits-1));
  return ans;
}

//input: a number n (to be converted to binary), the number of bits required to convert n to binary,
//and an address in memory to which to write the result
//conversion is done by lookup table
void getBinaryFromTable(unsigned long long n, unsigned short numBits, bool* ans)
{
  FILE* binTable;
  binTable = fopen("binTable.bin","rb");

  unsigned long long start = whereToLook(n);

  fseek(binTable,start,SEEK_SET);
  fread(ans,1,numBits,binTable);
  fclose(binTable);
}


//input: a number n (to be converted to binary), the number of bits required to convert n to binary,
//and an address in memory to which to write the result
void getBinary(unsigned long long n,unsigned short numBits,bool* binaryRep)
{

  cout << "GETTING NEW BINARY" << endl;
  if (numBits <= MAX_POWER)
  {
    getBinaryFromTable(n,numBits,binaryRep);
    return;
  }
  else
  {
    for (int i=0;i<numBits;i++)
      binaryRep[i]=false;

    unsigned short newNumBits = numBits;
    while (newNumBits > MAX_POWER)
    {
      cout << "...IN WHILE LOOP" << endl;
      n -= mypow(2,numBits-1);
      binaryRep[newNumBits-1] = true;
      newNumBits = getNumBitsNeeded(n);
    }
    getBinaryFromTable(n, newNumBits, binaryRep+(numBits-newNumBits));
    return;
  }


  for (int i=numBits - 1;i>=0;i--)
    if (n >= mypow(2,i))
    {
      binaryRep[i] = true;
      n -= mypow(2,i);
    }

  //  return binaryRep;
}


//generates a lookup table that contains the binary forms of all integers up to 2^MAX_POWER
//each binary number is stored with the minimum necessary disk space
void makeBinaryLookup()
{
  ofstream log("binlog.txt");
  FILE* binTable;
  binTable = fopen("binTable.bin", "wb");

  unsigned long long stop = mypow(2,MAX_POWER)-1;
  int percent=0;
  unsigned short buffer = 0;
  unsigned long long check = mypow(2,MAX_POWER-7);
  for (unsigned long long i = 0; i <= stop; i++)
  {
    if (i%check==0)
    {
      log << "done with " << percent << "/126" << endl;
      percent++;
    }

    buffer = getNumBitsNeeded(i);
    bool* arr;
    arr = new bool[buffer];
    getBinary(i,buffer,arr);
    fwrite(arr, 1, buffer, binTable);
    delete[] arr;
  }
  log.close();
  fclose(binTable);
  return;
  
}


//returns the negligibility of n
unsigned short getNeg(unsigned long long n)
{
  unsigned long long original = n;
  bool* binaryRep;

  int maxPower = 0;
  while (n >= pow(2,maxPower))
    maxPower++;

  //  cout<<"evaluated maxPower of " << original  << " to be " << maxPower << endl;

  binaryRep = new bool[maxPower];
  for (int i=0;i<maxPower;i++)
    binaryRep[i]=false;

  for (int i=maxPower - 1;i>=0;i--)
  {
    if (n >= pow(2,i))
    {
      binaryRep[i] = true;
      n -= pow(2,i);
    }
  }

  /*
  cout << "binary rep of " << original << " is: ";
  for (int i=maxPower-1;i>=0;i--)
    if (binaryRep[i])
      cout << "1";
    else
      cout << "0";
  cout << endl;
  */       

  unsigned short rank=0;
  unsigned short sum=0;
  for (int i=0;i<maxPower;i++)
    if (binaryRep[i])
    {
      rank+=1;
      sum+=i+1;
    }
  //  pair<int,int> ans(sum,rank);
  delete[] binaryRep;
  return rank+sum;
}

int main()
{
  // ******************* Generates the negligibility table:******************** 
  /*
  time_t start;
  double seconds;
  
  int NUM_POWERS = 23;
  int MAX_STACK = 15;

  ofstream log("log.txt");
  
  const unsigned long long size = pow(2,MAX_STACK);
  unsigned short a[size];

  FILE* negTable;
  negTable = fopen("negTableTest.bin", "wb");
  time(&start);

  unsigned long long buffsize = size*sizeof(unsigned short);
  unsigned long long stop = pow(2,NUM_POWERS-MAX_STACK);
  for (unsigned long long j = 0; j < stop; ++j)
  {
    unsigned long long stopplace = (j+1)*size;
    unsigned long long clumpsize = j*size;
    log << "Done with " << j / stop * 100 << "%  (" << j*size << " of " << stop*size << ")" << endl;
    //Some calculations to fill a[]
    for (unsigned long long i=j*size;i<stopplace;i++)
      a[i-clumpsize] = getNeg(i);

    fwrite(a, 1, buffsize, negTable);
  }
  time_t end;
  time(&end);
  seconds = difftime(end, start);
  cout << "time: " << seconds << endl;
  fclose(negTable);

  log << "TOOK " << seconds << " seconds" << endl;
  log.close();
  return 0;
*/
  // ******************* End negTable generation********************* 


  //  makeBinaryLookup();


  // Testing the conversion from decimal to binary...
  unsigned long long thing1 = mypow(2,14)+mypow(2,12)+mypow(2,10)+mypow(2,8);
  unsigned long long thing2 = mypow(2,50)+mypow(2,26)+mypow(2,10)+mypow(2,15);
  unsigned long long thing3 = 102859123841;
  unsigned short size1 = getNumBitsNeeded(thing1);
  unsigned short size2 = getNumBitsNeeded(thing2);
  unsigned short size3 = getNumBitsNeeded(thing3);
  bool* bin1; bin1 = new bool[size1]; getBinary(thing1,size1,bin1);
  bool* bin2; bin2 = new bool[size2]; getBinary(thing2,size2,bin2);
  for (int i=size1;i>=0;i--)
    if (bin1[i])
      cout << "1";
    else 
      cout << "0";
  cout << endl;
  for (int i=size2;i>=0;i--)
    if (bin2[i])
      cout << "1";
    else 
      cout << "0";
  cout << endl;


  bool* bin3; bin3 = new bool[size3]; getBinary(thing3,size3,bin3);
  
  for (int i=size3;i>=0;i--)
    if (bin3[i])
      cout << "1";
    else 
      cout << "0";
  cout << endl;

  /*
  cout << "mypow(2,30) = " << mypow(2,30) << endl;
  bool* bin = getBinary(mypow(2,30),MAX_POWER);
  cout << "binary of 2^30: ";
  for (int i=MAX_POWER;i>=0;i--)
    if (bin[i])
      cout << "1";
    else 
      cout << "0";
  cout << endl;
  */


  /*
  ofstream data("p_c.dat");

  int maxp_c = 30;
  unsigned long long* lessThan;
  unsigned long long* equalsTo;
  lessThan = getLessThanArray(maxp_c);
  equalsTo = getEqualsToArray(maxp_c);

  for (int i=0;i<maxp_c;i++)
    cout << i << " " << lessThan[i] << " " << equalsTo[i] << endl;

  data.close();
  delete[] lessThan;
  delete[] equalsTo;

  */
  //    cout << getNegFromTable(0);


  /*
  cout << "up to what p_c? ";
  int maxp_c;
  cin >> maxp_c;
  
  ofstream output("output.txt");

  pair<int,int>* negs;
  int ARRAY_SIZE = pow(2,20);
  negs = new pair<int,int>[ARRAY_SIZE];
  for (int i=0;i<ARRAY_SIZE;i++)
  {
    pair<int,int> cur = getNeg(i);
    negs[i] = cur;
  }
  
  for (int p_c = 2; p_c < maxp_c;p_c++)
  {
    int numLess = 0;
    int numEqual = 0;
    for (int i=1;i<pow(2,p_c-2)+1;i++)
    {
      pair<int,int> cur;
      if (i < pow(2,15))
	cur = negs[i];
      else
	cur = getNeg(i);
      
      if (cur.first + cur.second < p_c)
	numLess++;
      else if (cur.first + cur.second == p_c)
	numEqual++;
      
      //    double lv1 = lg2(i)+2;
      //double lv2 = lg2(i-pow(2,deg(i)))+deg(i)+4;
      //double lv3 = lg2(i - pow(2,deg(i)) - pow(2,deg(i-pow(2,deg(i))))) + deg(i-pow(2,deg(i))) + deg(i) + 6;
      //    double lv4 = lg2(i - pow(2,deg(i)) - pow(2,deg(i-pow(2,deg(i)))) - pow(2,deg(i-pow(2,deg(i))-pow(2,i-pow(2,i))))) 
      // + deg(i-pow(2,deg(i)) - pow(2,deg(i-pow(2,deg(i))))) + deg(i-pow(2,deg(i))) + deg(i) + 8;
      
      
      output << i << " " <<  p_c << " " << cur.first + cur.second << " " << numLess << " " << numEqual << endl;
      //	   << lv1 << " " 
      //	   << lv2 << " " 
      //	   << lv3 << " " 
      //	   << thing << endl;
      
    }
  }
  output.close();
  */
}

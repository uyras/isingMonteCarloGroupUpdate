#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

static default_random_engine generator;

static uint L = 100;
static int J = -1;
static uint n;
static int E;
static double T;

static unsigned long long int heatupMcSteps = 30000;
static unsigned long long int mcSteps = 400000;

static double sumSqr = 0, sum2Sqr = 0;

static vector<int> s, sMin;
static vector<uint> neigh;

inline int idx(const int i, const int j) { return i * L + j; }
inline int getI(const int idx) { return idx / L; }
inline int getJ(const int idx) { return idx % L; }

inline int PBC(const int i)
{
	if (i < 0) {
		return L + i;
	}
	else if (i >= L) {
		return i - L;
	}
	else {
		return i;
	}
}

int calcE() 
{
	E = 0;
    for (size_t i = 0; i < size_t(n); ++i) {
		E += -J * s[i] * s[neigh[i * 4 + 0]];
		E += -J * s[i] * s[neigh[i * 4 + 1]];
	}
	return E;
}

bool next() 
{
    for (size_t i = 0; i < size_t(n); ++i) {
		if (s[i] < 0) {
			s[i] *= -1;
			return true;
		}
		if (s[i] > 0)
			s[i] *= -1;
	}
	return false;
}

void mc(unsigned long long steps ){


    uniform_int_distribution<size_t> intrand(0, n - 1);
	uniform_real_distribution<double> prob(0, 1);


	int eOld = calcE();

    int eNew;
    size_t cand;
	double P1, P;

	for (unsigned long long step = 0; step < steps; ++step)
	{
		
		/*if (step % n == 0)
		{
			for (int i = 0; i < n - 1; ++i)
			{
				f << ((s[i] < 0) ? "-1," : "1,");
			}
			f << ((s[n - 1] < 0) ? "-1" : "1");
			f << endl;
		}*/

		cand = intrand(generator);

		eNew = eOld + J * 2 * s[cand] * (
			s[neigh[cand * 4 + 0]] +
			s[neigh[cand * 4 + 1]] +
			s[neigh[cand * 4 + 2]] +
			s[neigh[cand * 4 + 3]]
			);

		P = exp(-(eNew - eOld) / T);
		P1 = prob(generator);

        if (P1 < P) {
			s[cand] *= -1;
			eOld = eNew;
		}


		//сумма квадратов скользящего среднего значения Е
		sumSqr += (eOld - sumSqr)/double(step + 1);
		sum2Sqr += (eOld*eOld - sum2Sqr) / double(step +1);

	}
}

int main(int argc, char *argv[]) 
{
	if (argc!=3){
        cout<<"error reading console parameters."<<endl;
        cout<<"The format is as follows:"<<endl;
        cout<<argv[0]<<" <T> <rseed>"<<endl;
        cout<<"\t<T> - temperature, real number (>0 and <10)"<<endl;
        cout<<"\t<rseed> - random seed, integer (>1)"<<endl;
        return 0;
    }
    T = stod(argv[1]);
    uint rseed = stoi(argv[2]);

	generator.seed(rseed);

	n = L * L;
	
	s.resize(n, -1);

	sMin.resize(n, -1);
	neigh.resize(n * 4);

	for (int i = 0; i < L; ++i) 
    {
		for (int j = 0; j < L; ++j) 
		{
			int index = idx(i, j);
			neigh[index * 4 + 0] = idx(PBC(i - 1), PBC(j));
			neigh[index * 4 + 1] = idx(PBC(i), PBC(j + 1));
			neigh[index * 4 + 2] = idx(PBC(i + 1), PBC(j));
			neigh[index * 4 + 3] = idx(PBC(i), PBC(j - 1));
			
			if ((i % 2) ^ (j % 2))
			{
				s[index] = -1;
			}
			else
			{
				s[index] = +1;
			}
		}
	}

	sumSqr=0; sum2Sqr=0;
	mc(heatupMcSteps*n);

	sumSqr=0; sum2Sqr=0;
	mc(mcSteps*n);

	cout << T <<"\t" << sumSqr << "\t" << sum2Sqr << endl;

	return 0;
}

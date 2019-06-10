#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

static default_random_engine generator;

// setup variables
#define L 100
#define REPLICAS 10
static int J = -1;
static unsigned long long int heatupMcSteps = 3000; //number of steps to heatup
static unsigned long long int mcSteps = 4000; //number of steps to perform monte-carlo


// global variables
static uint N; //number of spins
static double eAvg = 0, e2Avg = 0; //averages
static double tMin, tMax; // temperatures for calculating B
static double B[REPLICAS]; // inverse temperatures for every replica
static int E[REPLICAS]; //local energy for every replica
static vector< vector<signed char> > s; //arrray of systems for each replica
static vector<uint> neigh; // array of neighbours
static vector<unsigned long long int> swapscount;

inline uint idx(const uint i, const uint j) { return i * L + j; }
inline uint getI(const uint idx) { return idx / L; }
inline uint getJ(const uint idx) { return idx % L; }

inline uint PBC(const int i)
{
    if (i < 0) {
        return L + uint(i);
    }
    else if (i >= int(L)) {
        return uint(i) - L;
    }
    else {
        return uint(i);
    }
}

int calcE(size_t replica)
{
    int E = 0;
    for (size_t i = 0; i < size_t(N); ++i) {
        E += -J * s[replica][i] * s[replica][neigh[i * 4 + 0]];
        E += -J * s[replica][i] * s[replica][neigh[i * 4 + 1]];
    }
    return E;
}

void mc(unsigned long long steps)
{
    int eTemp;
    size_t cand,replica2;
    double P1, P;

    uniform_int_distribution<size_t> intrand(0, N - 1);
    uniform_real_distribution<double> prob(0, 1);

    eTemp = calcE(0);
    for (size_t replica = 0; replica < REPLICAS; ++replica){
        E[replica] = eTemp; //same configs on all replicas
    }


    for (unsigned long long step = 0; step < steps; ++step)
    {
        for (size_t replica = 0; replica < REPLICAS; ++replica){

            cand = intrand(generator);

            eTemp = J * 2 * s[replica][cand] * (
                s[replica][neigh[cand * 4 + 0]] +
                s[replica][neigh[cand * 4 + 1]] +
                s[replica][neigh[cand * 4 + 2]] +
                s[replica][neigh[cand * 4 + 3]]
                );

            P = exp(-(eTemp) * B[replica]);
            P1 = prob(generator);

            if (P1 < P) {
                s[replica][cand] *= -1;
                E[replica] += eTemp;
            }


            if (replica==0){
                //сумма квадратов скользящего среднего значения Е
                eAvg += (E[replica] - eAvg)/double(step + 1);
                e2Avg += (E[replica]*E[replica] - e2Avg) / double(step +1);
            }

            if (prob(generator)<0.5){
                replica2=replica-1;
                if (replica2 > 0){
                    P = exp((B[replica2]-B[replica])*(E[replica2]-E[replica]));
                    P1 = prob(generator);
                    if (P1 < P){
                        s[replica].swap(s[replica2]);
                        swap(E[replica],E[replica2]);
                        ++swapscount[replica];
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc!=4){
        cout<<"error reading console parameters."<<endl;
        cout<<"The format is as follows:"<<endl;
        cout<<argv[0]<<" <T> <Tmax> <rseed>"<<endl;
        cout<<"\t<T> - temperature, real number (>0 and <10)"<<endl;
        cout<<"\t<Tmax> - maximal temperature, real number (>T and <15)"<<endl;
        cout<<"\t<rseed> - random seed, integer (>1)"<<endl;
        return 0;
    }

    // init seed
    uint rseed = uint(stoi(argv[3]));
    generator.seed(rseed);

    //init temperatures array
    tMin = stod(argv[1]);
    tMax = stod(argv[2]);

    N = L * L;

    for (size_t i=0; i<REPLICAS; ++i){ // inverse temperature as power function
        B[i] = 1./(pow(tMax/tMin,i/(double(REPLICAS-1)))*tMin);
    }

    s.resize(REPLICAS);
    for (size_t i=0; i<REPLICAS; ++i){ // inverse temperature as power function
        s[i].resize(N,-1);
    }

    neigh.resize(N * 4);

    for (uint i = 0; i < L; ++i)
    {
        for (uint j = 0; j < L; ++j)
        {
            uint index = idx(i, j);
            neigh[index * 4 + 0] = idx(PBC(int(i) - 1), PBC(int(j)));
            neigh[index * 4 + 1] = idx(PBC(int(i)), PBC(int(j) + 1));
            neigh[index * 4 + 2] = idx(PBC(int(i) + 1), PBC(int(j)));
            neigh[index * 4 + 3] = idx(PBC(int(i)), PBC(int(j) - 1));

            //set all replicas to ground states
            for (size_t i=0; i<REPLICAS; ++i){
                if ((i % 2) ^ (j % 2))
                    s[i][index] = -1;
                else
                    s[i][index] = +1;
            }
        }
    }

    swapscount.resize(REPLICAS,0);
    eAvg=0; e2Avg=0;
    mc(heatupMcSteps*N);
    //cout<<"swaps before resize: ";
    for (size_t i=0; i<REPLICAS; ++i){
        //cout<<swapscount[i]/double(heatupMcSteps*N)<<"; ";
        swapscount[i]=0;
    }
    //cout<<endl;

    eAvg=0; e2Avg=0;
    mc(mcSteps*N);
    /*cout<<"swaps after resize: ";
    for (size_t i=0; i<L; ++i){
        cout<<swapscount[i]/double(mcSteps*N)<<"; ";
    }
    cout<<endl;*/

    cout << tMin <<"\t" << eAvg << "\t" << e2Avg << endl;

    return 0;
}

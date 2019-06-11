#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

static default_random_engine generator;

// setup variables
#define L 100
#define REPLICAS 5
static int J = -1;
static unsigned long long int heatupMcSteps = 300; //number of steps to heatup
static unsigned long long int mcSteps = 400; //number of steps to perform monte-carlo


// global variables
static uint N; //number of spins
static double eAvg = 0, e2Avg = 0; //averages
static double tMin, tMax; // temperatures for calculating B
static double B[REPLICAS]; // inverse temperatures for every replica
static size_t BConnections[REPLICAS]; //connection between temperature and replica number
//for example, BConnections[1]=3 means that temperature B[1] is assigned to s[3] replica
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
    size_t cand,replica,replica2;
    double P1, P;

    uniform_int_distribution<size_t> intrand(0, N - 1);
    uniform_real_distribution<double> prob(0, 1);

    for (unsigned long long step = 0; step < steps; ++step)
    {
        // canonical MC step for every replica
        for (size_t BNum = 0; BNum < REPLICAS; ++BNum){

            replica = BConnections[BNum];

            cand = intrand(generator);

            eTemp = J * 2 * s[replica][cand] * (
                s[replica][neigh[cand * 4 + 0]] +
                s[replica][neigh[cand * 4 + 1]] +
                s[replica][neigh[cand * 4 + 2]] +
                s[replica][neigh[cand * 4 + 3]]
                );

            P = exp(-(eTemp) * B[BNum]);
            P1 = prob(generator);

            if (P1 < P) {
                s[replica][cand] *= -1;
                E[replica] += eTemp;
            }


            if (BNum==0){
                //сумма квадратов скользящего среднего значения Е
                eAvg += (E[replica] - eAvg)/double(step + 1);
                e2Avg += (E[replica]*E[replica] - e2Avg) / double(step +1);
            }
        }

        // PT exchange step
        for (size_t BNum = 1; BNum < REPLICAS; ++BNum){
            replica=BConnections[BNum-1];
            replica2 = BConnections[BNum];
            P = exp((B[BNum-1]-B[BNum])*(E[replica]-E[replica2]));
            P1 = prob(generator);
            if (P1 < P){
                swap(BConnections[BNum],BConnections[BNum-1]);
                ++swapscount[BNum];
            }
        }

        //check myself
        if (step==10000){
            for (replica = 0; replica < REPLICAS; ++replica){
                if (E[replica]!=calcE(replica))
                    cout<<replica<<": "<<E[replica]<<"!="<<calcE(replica)<<endl;
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
    tMin=0.001;
    tMax=4;

    N = L * L;

    /*for (size_t i=0; i<REPLICAS; ++i){ // inverse temperature as power function
        B[i] = 1./(pow(tMax/tMin,i/(double(REPLICAS-1)))*tMin);
        BConnections[i] = i;
    }*/

   double bb[] = {
       0.1       , 0.13233333, 0.22933333, 0.391     , 0.61733333,
              0.90833333, 1.264     , 1.68433333, 2.16933333, 2.719
    };
    for (size_t i=0; i<REPLICAS; ++i){ // inverse temperature as power function
           B[i] = 1./bb[i];
           BConnections[i] = i;
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
            //fill the array of neighbours for square lattice, with PBC
            uint index = idx(i, j);
            neigh[index * 4 + 0] = idx(PBC(int(i) - 1), PBC(int(j)));
            neigh[index * 4 + 1] = idx(PBC(int(i)), PBC(int(j) + 1));
            neigh[index * 4 + 2] = idx(PBC(int(i) + 1), PBC(int(j)));
            neigh[index * 4 + 3] = idx(PBC(int(i)), PBC(int(j) - 1));

            //set all replicas to ground states
            for (size_t k=0; k<REPLICAS; ++k){
                if ((i % 2) ^ (j % 2))
                    s[k][index] = -1;
                else
                    s[k][index] = +1;
            }
        }
    }

    //calc initial energies for replicas
    int eTemp = calcE(0);
    for (size_t replica = 0; replica < REPLICAS; ++replica){
        E[replica] = eTemp; //same configs on all replicas
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

    cout << 1./B[0] <<"\t" << eAvg << "\t" << e2Avg << "\t" << rseed << endl;

    //check myself again
    cout<<"# Acceptance rates:"<<endl;
    for (size_t i=0; i<REPLICAS; ++i){
        cout<<"# "<<1./B[i]<<"J:"<< swapscount[i]/double(mcSteps*N) <<"; "<<endl;;
    }

    return 0;
}

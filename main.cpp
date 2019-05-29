#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <fstream>
#include <random>
#include <bitset>

using namespace std;
static default_random_engine generator;

#define L 100
#define N L*L
#define NEIGHBOURS 4
#define KERNEL_L 3 // linear size of core
#define KERNEL_N KERNEL_L*KERNEL_L // number of spins in core
#define BORDER_N (KERNEL_L*4) //number of neighbouring surround around of core
#define BORDER_L (KERNEL_L+2) //linear size of border
#define BLOCK_N (BORDER_L*BORDER_L) // total spinsin kernel + in border, plus 4 unused spins in edge
#define J 1

#define nemax (4*N)

static int rseed = 1;
static unsigned long long int heatupMcSteps = 3000;
static unsigned long long int mcSteps = 400000;
static double t=0.001;
static int energy;
static double eAverage=0,e2Average=0;


struct EMtype {
    int E;
    int M;
    vector <uint> states;

    EMtype(int e, int m): E(e), M(m){}
    EMtype(int e, int m, uint s): E(e), M(m){states.push_back(s);}
    inline bool operator==(const EMtype &a) const {return E==a.E&&M==a.M;}
};

struct Etype {
    int E;
    vector <uint> states;

    Etype(int e): E(e){}
    Etype(int e, uint s): E(e) {states.push_back(s);}
    inline bool operator==(const Etype &a) const {return E==a.E;}
};

struct pairtype {
    int ai;
    int aj;
    int bi;
    int bj;
};

struct postype {
    int i,j;
};

struct simpleEM {
    int E;
    int M;
    double prob;
};

struct averageValType {
    double avgLocE;
    double avgLocE2;
    double totalProbability;
    int minimalEnergy;
};

static signed char s[N]; // spin values over all system
static int neigh[N*NEIGHBOURS];

static vector< vector <simpleEM> > states; //access through states[b][i].{E|M|prob}
// where b is bit-typed configuration of border, i - sequential configuration of kernel;

static vector< averageValType > coreAverages; // array of minimal energies of kernel grouped by border.

static vector<pairtype> pairmask; // mask for all pair interactions counting from top left corner of core spin
static postype borderpos[BORDER_N];



const int LSHIFT = L*10; //shift of L along the numbers line to get pbc() function working

inline int idx(int _i, int _j){return _i*L+_j;}
inline int geti(int _n){return _n/L;}
inline int getj(int _n){return _n%L;}
//inline int pbc(int _n){return (_n<0)?(L*100+_n)%L:(_n>=L)?_n%L:_n;} // replaced for faster code
inline int pbc(int _n){return (L*LSHIFT+_n)%L;}


inline int block_idx(int _i, int _j){return _i*BORDER_L+_j;}
inline int block_geti(int _n){return _n/BORDER_L;}
inline int block_getj(int _n){return _n%BORDER_L;}

void mc(unsigned long long int );
void step();
void single();
void showLocSys(signed char *locs, size_t LL);

void spinset(){
    int index = 0;
    for (int i=0; i<L; ++i){
        for (int j=0; j<L; ++j){
            index = idx(i,j);
            //s[index] = (lrand48()%2==0)?-1:+1;
            s[index] = 1;

            neigh[index*NEIGHBOURS + 0] = idx(i,pbc(j+1)); //right
            neigh[index*NEIGHBOURS + 1] = idx(pbc(i+1),j); //bottom
            neigh[index*NEIGHBOURS + 2] = idx(i,pbc(j-1)); //left
            neigh[index*NEIGHBOURS + 3] = idx(pbc(i-1),j); //top
        }
    }
}

void genmask(){
    //gen upper pairs
    for (int i=0; i<KERNEL_L; ++i){
        pairmask.push_back({0,i,-1,i}); //upper pairs
        pairmask.push_back({KERNEL_L-1,i,KERNEL_L,i}); //bottom pairs
        pairmask.push_back({i,0,i,-1}); //left pairs
        pairmask.push_back({i,KERNEL_L-1,i,KERNEL_L}); //right pairs
    }
    for (int i=0; i<KERNEL_L; ++i){
        for (int j=0; j<KERNEL_L; ++j){
            if (j<KERNEL_L-1)
                pairmask.push_back({i,j,i,j+1}); //right pairs
            if (i<KERNEL_L-1)
                pairmask.push_back({i,j,i+1,j}); //bottom pairs
        }
    }

    int p=0;
    //fill the right edge
    for (int i=0; i<KERNEL_L; ++i){
        borderpos[p] = {i,KERNEL_L}; ++p;
    }

    //fill the bottom edge
    for (int i=0; i<KERNEL_L; ++i){
        borderpos[p] = {KERNEL_L,i}; ++p;
    }

    //fill the left edge
    for (int i=0; i<KERNEL_L; ++i){
        borderpos[p] = {i,-1}; ++p;
    }

    //then fill the top edge
    for (int i=0; i<KERNEL_L; ++i){
        borderpos[p] = {-1,i}; ++p;
    }
}

void gatherConf(){
    const int borderConfigMax=1<<BORDER_N, kernelConfigMax=1<<KERNEL_N;

    states.resize(borderConfigMax);
    coreAverages.resize(borderConfigMax,{0,0,0,99999});

    signed char locs[BLOCK_N]; // local lattice of spins, with borders

    int index = 0;
    uint tempBorderConfigNum=0, tempKernelConfigNum=0, kernelConfigNum=0;
    int totalM, tempM, tempE;

    //enumerate surroundings
    for (uint borderConfigNum=0; borderConfigNum<borderConfigMax; ++borderConfigNum ){
        states[borderConfigNum].resize(kernelConfigMax);
        totalM=0;
        //fill config bits in array
        tempBorderConfigNum=borderConfigNum;

        //fill edges
        for (int i=0; i<BORDER_N; ++i){
            index=block_idx(borderpos[i].i+1,borderpos[i].j+1);
            totalM += (locs[index]=(tempBorderConfigNum&1)?1:-1);
            tempBorderConfigNum = tempBorderConfigNum>>1;
        }

        // enumerate configs
        for ( kernelConfigNum=0; kernelConfigNum < kernelConfigMax; ++kernelConfigNum ){
            tempM=totalM;
            tempKernelConfigNum=kernelConfigNum;
            //fill the center square
            for (int i=1; i<KERNEL_L+1; ++i){
                for (int j=1; j<KERNEL_L+1;++j){
                    index = block_idx(i,j);
                    tempM += (locs[index] = (tempKernelConfigNum&1)?1:-1);
                    tempKernelConfigNum = tempKernelConfigNum>>1;
                }
            }

            //get the energy
            tempE=0;
            for (ulong i=0; i<pairmask.size(); ++i){
                tempE +=
                        locs[block_idx(1+pairmask[i].ai,1+pairmask[i].aj)] *
                        locs[block_idx(1+pairmask[i].bi,1+pairmask[i].bj)];
            }
            tempE *= -J;

            //add value to dataset
            states[borderConfigNum][kernelConfigNum] = {tempE,tempM,0};

            // store the minimal energy
            if (tempE < coreAverages[borderConfigNum].minimalEnergy)
                coreAverages[borderConfigNum].minimalEnergy = tempE;
        }

        // calc thermal averages
        for ( kernelConfigNum=0; kernelConfigNum < kernelConfigMax; ++kernelConfigNum ){
            states[borderConfigNum][kernelConfigNum].prob =
                    exp(-(states[borderConfigNum][kernelConfigNum].E-coreAverages[borderConfigNum].minimalEnergy)/t);

            coreAverages[borderConfigNum].totalProbability += states[borderConfigNum][kernelConfigNum].prob;

            coreAverages[borderConfigNum].avgLocE += states[borderConfigNum][kernelConfigNum].E*
                                                     states[borderConfigNum][kernelConfigNum].prob;

            coreAverages[borderConfigNum].avgLocE2 += states[borderConfigNum][kernelConfigNum].E*
                                                      states[borderConfigNum][kernelConfigNum].E*
                                                      states[borderConfigNum][kernelConfigNum].prob;
        }
        coreAverages[borderConfigNum].avgLocE /= coreAverages[borderConfigNum].totalProbability;
        coreAverages[borderConfigNum].avgLocE2 /= coreAverages[borderConfigNum].totalProbability;
    }
}

void energyset(){
    energy = 0;
    for (int i=0; i<N; ++i){
        for (int j=0; j<NEIGHBOURS; ++j){
            energy += s[i]*s[neigh[i*NEIGHBOURS+j]];
        }
    }
    energy *= J*-0.5;
}

uint getBorderConf(const int di, const int dj){
    uint tempC=0, tbit;
    for (int i=BORDER_N-1; i>=0; --i){
        tbit = (s[idx(pbc(borderpos[i].i+di),pbc(borderpos[i].j+dj))]==-1)?0:1;
        if (i>0)
            tempC = (tempC | tbit)<<1;
        else
            tempC = tempC | tbit;
    }
    return tempC;
}

uint getKernelConf(const int di, const int dj){
    uint tempC=0, tbit;
    for (int i=KERNEL_L-1; i>=0; --i){
        for (int j=KERNEL_L-1; j>=0;--j){
            tbit = (s[idx(pbc(i+di),pbc(j+dj))]==-1)?0:1;
            if (i>0 || j>0)
                tempC = (tempC | tbit)<<1;
            else
                tempC = tempC | tbit;
        }
    }
    return tempC;
}

void setKernelConf(uint c, const int di, const int dj){
    //fill the center square
    for (int i=0; i<KERNEL_L; ++i){
        for (int j=0; j<KERNEL_L;++j){
            s[idx(pbc(i+di),pbc(j+dj))] = (c&1)?1:-1;
            c=c>>1;
        }
    }
}

void showSys(){
    cout<<"=============================="<<endl;
    for (int i=0; i<L; ++i){
        for (int j=0; j<L; ++j){
            cout<<((s[idx(i,j)]>0)?"+":"-");
        }
        cout<<endl;
    }
    cout<<"=============================="<<endl;
}

void showLocSys(signed char *locs, size_t LL){
    cout<<"=============================="<<endl;
    for (size_t i=0; i<LL; ++i){
        for (size_t j=0; j<LL; ++j){
            cout<<((locs[i*LL+j]>0)?"+":"-");
        }
        cout<<endl;
    }
    cout<<"=============================="<<endl;
}

void mc(unsigned long long int steps=1)
/*
        monte carlo update
*/
{
    unsigned long long int n;

    uint totalKernelConfigs=1<<KERNEL_N;

    int coreI, coreJ;
    int coreSpin;
    uniform_int_distribution<int> selectSpinRandom(0,N-1);

    uint kernelConfigNum;

    for( n = 0; n <= steps; ++n){


        int eShift=0;

        coreSpin = selectSpinRandom(generator);
        coreI = geti(coreSpin);
        coreJ = getj(coreSpin);

        // get configurations of kernel and border
        uint oldKernelConfig=getKernelConf(coreI,coreJ);
        uint borderConfig = getBorderConf(coreI,coreJ);


        eShift = energy - states[borderConfig][oldKernelConfig].E;


        /*energyset();
        if (eShift+currE->E != energy){
            cout<<"E error again~!"<<endl;
        }*/

        //строим локальную статсумму
        double eAverageLocGlob=0, e2AverageLocGlob=0;//попутно считаем средние параметры

        eAverageLocGlob  = eShift + coreAverages[borderConfig].avgLocE;
        e2AverageLocGlob = eShift*eShift + 2*eShift*coreAverages[borderConfig].avgLocE + coreAverages[borderConfig].avgLocE2;


        // Update total averages (sliding average value)
        eAverage  += (eAverageLocGlob  - eAverage )/(n+1);
        e2Average += (e2AverageLocGlob - e2Average)/(n+1);



        // выбираем нового кандидата
        uniform_real_distribution<double> realDistribution(0,coreAverages[borderConfig].totalProbability);
        double rval = realDistribution(generator);
        double tval=0;
        for (kernelConfigNum=0; kernelConfigNum<totalKernelConfigs; ++kernelConfigNum){
            if (rval>=tval && rval<tval+states[borderConfig][kernelConfigNum].prob){
                break;
            } else {
                tval+=states[borderConfig][kernelConfigNum].prob;
            }
        }
        setKernelConf(kernelConfigNum,coreI,coreJ);
        energy = eShift+states[borderConfig][kernelConfigNum].E;
    }
}

int main(int argc, char *argv[])
{
    //cout<<mcSteps<<endl;
    if (argc!=3){
        cout<<"error reading console parameters."<<endl;
        cout<<"The format is as follows:"<<endl;
        cout<<argv[0]<<" <T> <rseed>"<<endl;
        cout<<"\t<T> - temperature, real number (>0 and <10)"<<endl;
        cout<<"\t<rseed> - random seed, integer (>1)"<<endl;
        return 0;
    }
    t = stod(argv[1]);
    rseed = stoi(argv[2]);


    generator.seed(rseed);

    //cout<<"started"<<endl;
    spinset();
    genmask();
    gatherConf();
    energyset();

    //cout<<"mask created"<<endl;
    //cout<<energy<<endl;
    mc(heatupMcSteps*N);
    eAverage=0;
    e2Average=0;
    mc(mcSteps*N);

    cout<<t<<"\t"<<eAverage<<"\t"<<e2Average<<"\t"<<rseed<<endl;

    return 0;

    /*showSys();

    for (int n=0; n<100; ++n){

        //get new config of kernel
        uint newconf=lrand48() % (1<<KERNEL_N), oldconf;

        int pi=lrand48() % L, pj=lrand48() % L;
        oldconf = getKernelConf(pi,pj);
        int oldE=-9999;

        //get edges
        uint tempC=getEdgesConf(pi,pj);
        //cout<<tempC<<endl;

        for (auto it2: states[tempC]){
            for (auto it3: it2.states){
                if (it3==oldconf){
                    oldE = it2.E;
                    break;
                }
            }
        }

        cout<<"n="<<n<<", pi="<<pi<<", pj="<<pj<<", edges="<<bitset<BORDER_N>(tempC)<<", oldConf="<<bitset<KERNEL_N>(oldconf)<<", newConf="<<bitset<KERNEL_N>(newconf)<<endl;

        for (auto it2: states[tempC]){
            for (auto it3: it2.states){
                if (it3==newconf){
                    cout<<"\t"<<energy<<"\t"<<it2.E<<"\t"<<oldE<<"\t"<<energy+it2.E-oldE<<"\t";
                    break;
                }
            }
        }


        setKernelConf(newconf,pi,pj);
        energyset();
        cout<<energy<<endl;
        showSys();
    }


    for (ulong kk=0; kk<(1<<BORDER_N); ++kk){
        cout<<kk<<" ("<<bitset<BORDER_N>(kk)<<") config:"<<endl;
        for (auto it = states[kk].begin(); it!=states[kk].end(); ++it) {
            cout<<it->E<<":\t";
            for (auto it2: it->states){
                cout<<it2<<"("<<bitset<KERNEL_N>(it2)<<"),";
            }
            std::cout << endl;
        }
        cout<<"---------------------------------"<<endl;
    }

    for (ulong kk1=0; kk1<(1<<BORDER_N); ++kk1){
        cout<<kk1<<" ("<<bitset<BORDER_N>(kk1)<<") pairs: ";
        for (ulong kk2=0; kk2<(1<<BORDER_N); ++kk2){
            if (kk1==kk2) continue;
            bool cmpflag=true;
            for (ulong it1 = 0; it1 < states[kk1].size(); ++it1) {
                if (states[kk1][it1].E == states[kk2][it1].E && cmpflag){
                    for (ulong it2 =0; it2 < states[kk1][it1].states.size(); ++it2){
                        if (states[kk1][it1].states[it2]!=states[kk2][it1].states[it2]){
                            cmpflag=false;
                            break;
                        }
                    }
                } else {
                    cmpflag=false;
                    break;
                }
            }

            if (cmpflag)
                cout<<kk2<<"("<<bitset<BORDER_N>(kk2)<<"), ";
        }
        cout<<endl;
    }

    return 0;*/
}

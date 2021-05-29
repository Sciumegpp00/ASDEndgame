#define DEBUG

#include <iostream>
#include <stdlib.h>
//#include "lib/endgame.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>

#define INPUT_FILE "../inputs/input1.txt"
#define OUTPUT_FILE "outputs/output0.txt"

using namespace std;

struct Rock {
    int mass;
    int energy;

    void printRock(ofstream *out) {
        *out << mass << " " << energy << endl;
    }
};

struct City {
    vector<Rock> rocks;
};

//struct Edge {
//    int weight;
//    short c1;
//    short c2;
//};

int getWeight(short city1, short city2);
void input();
void output();
void bbTsp(vector<short>* path, int cost, vector<short>* remaining, int n, int i);
void removeVector(vector<short>* v, short x);


short nCitiesTot;
short nDifferentRocks;
int nTakenRocks; //FIXME non dovrebbe servire perché tanto dobbiamo stampare tutti i numeri delle pietre anche quelle non prese con -1
short source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
int **weights;
int globalMinWeight;
int minCost; //costo di un percorso minimo a caso (il primo percorso che prendo)
vector<short>* minSol = nullptr;
vector<short>* remainingBbTsp = nullptr;
//int edgeCounts;
//Edge* edges;

//int edgeCmp(Edge x, Edge y) {
//    return x.weight < y.weight
//           ? -1
//           : x.weight > y.weight;
//}

int main() {
    input();
    /*
    auto err = heapsort(edges, sizeof(Edge), edgeCounts, reinterpret_cast<int (*)(const void *, const void *)>(&edgeCmp));
    if(err == -1) {
        cout << "Error in heapsort!\n";
    }*/

    vector<short> sol;
    sol.push_back(source);

    bbTsp(&sol, 0, remainingBbTsp, nCitiesTot-1, 1);

    for(int i = 0; i < minSol->size(); i++){
        cout << (*minSol)[i] << " ";
    }
    cout << source << endl;
    cout << "Cost: " << minCost << endl;

//    for (u_long i = minSol->size()-1; i >= 0; i--) {
//    }

    //output();
    return 0;
}

void input() {
    ifstream in(INPUT_FILE);

    in >> nCitiesTot >> source;
    in >> nDifferentRocks >> C >> R >> vmin >> vmax;
    //cout << nCitiesTot << " " << source << endl << nDifferentRocks << " " << C << " " << R << " " << vmin << " " << vmax << endl;

    rocks = new Rock[nDifferentRocks];

    for (int i = 0; i < nDifferentRocks; ++i) {
        in >> rocks[i].mass >> rocks[i].energy; //mass and energy for each rock
        //rocks[i].printRock(&cout);
    }

    int rocksPerCity;
    int city;
    cities = new City[nCitiesTot];

    for (int i = 0; i < nDifferentRocks; ++i) { //i-th rock
        in >> rocksPerCity;
        //cout << rocksPerCity << endl;

        for (int j = 0; j < rocksPerCity; ++j) {
            in >> city;

            cities[city].rocks.push_back(rocks[i]);
            //cout << city << " ";
        }
    }

    weights = new int *[nCitiesTot - 1]; //1 is 0
    remainingBbTsp = new vector<short>();

//    edgeCounts = nCitiesTot * (nCitiesTot - 1) / 2;
//    edges = new Edge[edgeCounts];
    globalMinWeight = INT_MAX;

    //inserimento dei pesi tra città
    for(short i=0; i<nCitiesTot-1; i++){
        weights[i] = new int[i+1];
        for(short j=0; j<i+1; j++){
            in >> weights[i][j];

            if(weights[i][j] < globalMinWeight)
                globalMinWeight = weights[i][j];
        }
    }

    minSol = new vector<short>();
    minSol->push_back(source);
    minCost = 0;

    for (short i = 0; i < nCitiesTot; ++i) {
        if(i != source) { //FIXME ma se va da 0 a n-2 manca una città
            remainingBbTsp->push_back(i);

            // Creating a minimum path (random)
            minCost += getWeight(i, minSol->back()); //sum the cost to the minimum
            minSol->push_back(i); //add the new node in the solution
        }
    }
    minCost += getWeight(minSol->back(), source);
}

// TODO fix total variables
void output() {
    ofstream out(OUTPUT_FILE);

    //final glove's energy, rocks' energy and time for the journey
//    out << scientific << setprecision(10) << E << " "; //total glove's energy
//    out << scientific << setprecision(10) << G << " "; //total rocks energy
//    out << scientific << setprecision(10) << T << endl; //time used

    //print all taken rocks
    for (int i = 0; i < nTakenRocks; i++) { //FIXME n totali di pietre --> se non viene presa una viene stampato -1
        printf("%d ", rocks[i]);
    }
    printf("\n");

    //print the path (ordered, t0 = tN = source)
    for(int i = 0; i < minSol->size(); i++){
        cout << (*minSol)[i] << " ";
    }
    cout << source << endl;

    //TODO questo sotto non serve più giusto?????
//    for (int i = 0; i < nCitiesTot; i++) {
//        cities[i].printCity(&out);
//    }
    printf("\n***\n");
}

int calcLb(short origin, vector<short>* remaining, int cost) {
    int out, back, transfer, costLocal, minTransferOut, localTransferOut;

    out = INT_MAX;
    back = INT_MAX;
//    transfer = INT_MAX;
//    out = getWeight(origin, (*remaining)[0]);
//    back = getWeight((*remaining)[0], source);
//    transfer = out + getWeight((*remainingBbTsp)[0], (*remainingBbTsp)[1]);

    for (int i = 0; i < remaining->size(); i++) {
        // calculate minimum cost of exiting edges at this point
        costLocal = getWeight(origin, (*remaining)[i]);
        if(costLocal < out) {
            out = costLocal;
        }

        // Calculate for each city the best way to pass into it (best duo entry + exit)
        // origin --> qualsiasi nodo + qualsiasi nodo
//        minTransferOut = INT_MAX;
//        for (int j = 0; j < remaining->size(); j++){ //consider the remainingBbTsp nodes
//            if(i != j){
//                localTransferOut = getWeight((*remaining)[i], (*remaining)[j]);
//                if(localTransferOut < minTransferOut)
//                    minTransferOut = localTransferOut;
//            }
//        }
//        transfer += minTransferOut;

        //calculate lb cost to go back to source of the nearest cities
        costLocal = getWeight((*remaining)[i], source);
        if(costLocal < back) {
            back = costLocal;
        }
    }

    transfer = (globalMinWeight * remaining->size());
    // TODO: Try to remove /2 for remainingBbTsp nodes > 10 or 15 for example
//    return cost + ((out + back + transfer) / 2);
    return cost + out + back + transfer;
}

void bbTsp(vector<short>* path, int cost, vector<short>* remaining, int n, int i) {
    short origin = path->back();
    vector<short> choices(*remaining);
    int lb, localCost;

    for (auto it = choices.begin() ; it != choices.end(); ++it) {
        localCost = cost + getWeight(origin, *it);
        path->push_back(*it);
        removeVector(remaining, *it);

        if(i < n) {
            lb = calcLb(*it, remaining, localCost);
            if(lb < minCost){
                bbTsp(path, localCost, remaining, n, i + 1);
            }
        } else {
            cost = localCost;
            cost += getWeight(*it, (*path)[0]); //path->front() Returns a reference to the first element in the vector.
            if(cost < minCost) {
                if(minSol) {
                    delete minSol;
                    minSol = nullptr;
                }
                minSol = new vector<short>(*path);
                minCost = cost;
            }
        }

        path->pop_back();
        remaining->push_back(*it);
    }
}

int getWeight(short city1, short city2) { //city 2 is the line and city 1 is the arrival city (bigger)
    return city1 < city2
           ? weights[city2 - 1][city1]
           : weights[city1 - 1][city2];
}

int calcV(int vMax, int vMin, int W, int C) {
    return vmax - W * ((vMax - vMin) / C);
}

int calcEnergy(int gloveEnergy, int totalTime){
    return gloveEnergy - R * totalTime;
}

void removeVector(vector<short>* v, short x) {
    for (auto it = v->begin() ; it != v->end(); it++) {
        if(*it == x) {
            v->erase(it);
            return;
        }
    }
}

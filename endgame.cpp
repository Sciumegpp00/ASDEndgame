#define DEBUG

#include <iostream>
#include <stdlib.h>
#include "lib/endgame.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>

#define INPUT_FILE "inputs/input0.txt"
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
int nTakenRocks;
short source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
int **weights;
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

    bbTsp(minSol, 0, remainingBbTsp, nCitiesTot - 1, 1);

    

    output();
    return 0;
}

void input() {
    ifstream in(INPUT_FILE);

    in >> nCitiesTot >> source;
    in >> nDifferentRocks >> C >> R >> vmin >> vmax;
    //out << nCitiesTot << " " << source << endl << nDifferentRocks << " " << C << " " << R << " " << vmin << " " << vmax << endl;

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
    remainingBbTsp = new vector<short>(nCitiesTot - 1);

//    edgeCounts = nCitiesTot * (nCitiesTot - 1) / 2;
//    edges = new Edge[edgeCounts];
    minSol = new vector<short>(nCitiesTot);
    minSol->push_back(source);
    minCost = 0;

    for (short i = 1; i < nCitiesTot; ++i) {
        if(i != source) {
            remainingBbTsp->push_back(i);
        }

        weights[i] = new int[i];
        for (short j = 0; j < i; ++j) {
            in >> weights[i][j];

            //creating a minimum path (random)
            if(j == minSol->back()){ //if we're considering the distance between the previous point and another (anyone)
                minCost += weights[i][j]; //sum the cost to the minimum
                minSol->push_back(i); //add the new node in the solution
            }
//            edges[(i * (i - 1) / 2) + j].weight = weights[i][j];
//            edges[(i * (i - 1) / 2) + j].c1 = i;
//            edges[(i * (i - 1) / 2) + j].c2 = j;
            //cout << weights[i][j] << " ";
        }
        //cout << endl;
    }
}

// TODO fix total variables
/*void output() {
    ofstream out(OUTPUT_FILE);

    //final glove's energy, rocks' energy and time for the journey
    out << scientific << setprecision(10) << E << " "; //total glove's energy
    out << scientific << setprecision(10) << G << " "; //total rocks enrgy
    out << scientific << setprecision(10) << T << endl; //time used

    //print all taken rocks
    for (int i = 0; i < nTakenRocks; i++) {
        printf("%d ", rocks[i]);
    }
    printf("\n");

    //print the path (ordered, t0 = tN = source)
    for (int i = 0; i < nCitiesTot; i++) {
        cities[i].printCity(&out);
    }
    printf("\n***\n");
}*/

int calcLb(short origin, vector<short>* remaining, int cost) {
    int out, back, transfer, costLocal;

    out = getWeight(origin, (*remaining)[0]);
    back = getWeight((*remaining)[0], source);
//    transfer = out + getWeight((*remainingBbTsp)[0], (*remainingBbTsp)[1]);
    for (int i = 1; i < remaining->size(); i++) {
        // calculate minimum cost of exiting edges at this point
        costLocal = getWeight(origin, (*remaining)[i]);
        if(costLocal < out) {
            out = costLocal;
        }

        // Calculate for each city the best way to pass into it (best duo entry + exit)
        // origin --> qualsiasi nodo + qualsiasi nodo
//        int bestExit = weights[(*remainingBbTsp)[i]][(*remainingBbTsp)[i-1]];
//        for (int j = 0; j < remainingBbTsp->size(); j++){ //consider the remainingBbTsp nodes
//            if(i != j){
//
//            }
//        }

        //calculate lb cost to go back to source of the nearest cities
        costLocal = getWeight((*remaining)[i], source);
        if(costLocal < back) {
            back = costLocal;
        }
    }

    // TODO: Try to remove /2 for remainingBbTsp nodes > 10 or 15 for example
    return cost + out + back;
}


void bbTsp(vector<short>* path, int cost, vector<short>* remaining, int n, int i) {
    short origin = path->back();
    vector<short> choices = *remaining;
    int lb, arriveCost;

    for (auto it = choices.begin() ; it != choices.end(); ++it) {
        path->push_back(*it);
        removeVector(remaining, *it);

        if(i < n) {
            arriveCost = getWeight(origin, *it);
            lb = calcLb(*it, &choices, cost + arriveCost);
            if(lb < minCost){
                bbTsp(path, cost + arriveCost, remaining, n, i + 1);
            }
        } else {
            cost += getWeight(*it, (*path)[0]); //path->front() Returns a reference to the first element in the vector.
            if(minSol) {
                delete minSol;
                minSol = nullptr;
            }
            minSol = new vector<short>(choices);
            minCost = lb;
        }

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
    for (auto it = v->begin() ; it != v->end(); ++it) {
        if(*it == x)
            v->erase(it);
    }
}

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

// @tocompile /usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -static -s -o endgame grader.cpp endgame.cpp
// problem with compiling with Clion --> modify cMakeList.txt
struct Rock {
    int mass;
    int energy;

    void printRock(ofstream *out) {
        *out << mass << " " << energy << endl;
    }
};

struct City {
    int n;
    vector<Rock> rocks;

    void printCity(ofstream *out) {
        *out << n << endl;
    }
};

//struct Edge {
//    int weight;
//    short c1;
//    short c2;
//};


short int nCitiesTot;
short int nDifferentRocks;
int nTakenRocks;
int source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
int **weights;
//int edgeCounts;
//Edge* edges;

int getWeight(short city1, short city2);
void input();
void output();
void bbTps();

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
    bbTps();

    output();
    return 0;
}

void removeVector(vector<short>* v, short x) {
    for (auto it = v->begin() ; it != v->end(); ++it) {
        if(*it == x)
            v->erase(it);
    }
}

int calcLb(short origin, vector<short>* choices, int cost) {
    int lb, out, back, *transfer, costLocal;
    short outCity;
    double minLb;

    if(choices->empty())
        return 0; // FIXME: Add cost to come back

    //calculate minimum cost of exiting edges at this point
    out = getWeight(origin, (*choices)[0]);
    outCity = 0;
    for (int i = 1; i < choices->size(); i++) {
        costLocal = getWeight(origin, (*choices)[i]);
        if(costLocal < out) {
            out = costLocal;
            outCity = (*choices)[i];
        }
    }

    //calculate lb cost to go back to source of the nearest cities
    back = getWeight(origin, (*choices)[0]);
    for (int i = 1; i < choices->size(); i++) {
        costLocal = getWeight((*choices)[i], (*choices)[0]);
        if(costLocal < back) {
            back = costLocal;
        }
    }

    //calculate min future path
    minLb = cost + ((out+back + transfer[0])/2);
    for (int i = 1; i < nCitiesTot; i++) {
        lb = cost + ((out+back + transfer[i])/2);
        if(lb < minLb) minLb = lb;
    }

    return minLb;
}

void bbTsp(vector<short>* path, int cost, vector<short>* leftChoices, int n, int i) {
    vector<short> choices = *leftChoices;
    vector<short> minSol;
    int minCost = 100000000; //FIXME: initial cost
    int lb;

    for (auto it = choices.begin() ; it != choices.end(); ++it){
        path->push_back(*it);
        removeVector(leftChoices, *it);

        if(i < n){
            lb = calcLb(*it, &choices, cost);
            if(lb < minCost){
                bbTsp(path, cost + getWeight((short) *(it-1), (short) *it), leftChoices, n, i+1);
            }
        }else{
            cost += getWeight((short) *(it), (short) (*path)[0]); //path->front() Returns a reference to the first element in the vector.
            vector<short> minSol = choices;
            minCost = lb;
        }
    }
}

int getWeight(short city1, short city2) { //city 2 is the line and city 1 is the arrival city (bigger)
    return city1 < city2
           ? weights[city2 - 1][city1]
           : weights[city1 - 1][city2];
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

            cities[city].rocks.push_back(rocks[i]); //push(number of rock)
            //cout << city << " ";
        }
    }

    weights = new int *[nCitiesTot - 1]; //1 is 0

//    edgeCounts = nCitiesTot * (nCitiesTot - 1) / 2;
//    edges = new Edge[edgeCounts];

    for (short i = 1; i < nCitiesTot; ++i) {
        weights[i] = new int[i];
        for (short j = 0; j < i; ++j) {
            in >> weights[i][j];
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

int calcV(int vMax, int vMin, int W, int C) {
    return vmax - W * ((vMax - vMin) / C);
}

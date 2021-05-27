#include <iostream>
#include <stdlib.h>
#include "endgame.h"
#include <iomanip>
#include <fstream>
#include <vector>

// @tocompile /usr/bin/g++ -DEVAL -std=c++11 -O2 -pipe -static -s -o endgame grader.cpp endgame.cpp
// problem with compiling with Clion --> modify cMakeList.txt
struct Rock {
    int mass;
    int energy;

    void printRock(ofstream* out){
        *out << mass << " " <<  energy << endl;
    }
};

struct City {
    int n;
    vector<Rock> rocks;

    void printCity(ofstream* out){
        *out << n << endl;
        for (int i = 0; i < rocks.size(); ++i) {
            *out << rocks[i].mass << " " <<  rocks[i].energy;
        }
        *out << endl;
    }
};


int nCitiesTot;
int nDifferentRocks;
int nTakenRocks;
int source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
int **weights;

int getWeight(int city1, int city2) { //city 2 is the line and city 1 is the arrival city (bigger)
    return city1 < city2
           ? weights[city2 - 1][city1]
           : weights[city1 - 1][city2];
}

void input() {
    ifstream in("input0.txt");
    ofstream out("output.txt");

    in >> nCitiesTot >> source;
    in >> nDifferentRocks >> C >> R >> vmin >> vmax;
    //out << nCitiesTot << " " << source << endl << nDifferentRocks << " " << C << " " << R << " " << vmin << " " << vmax << endl;

    rocks = new Rock[nDifferentRocks];

    for (int i = 0; i < nDifferentRocks; ++i) {
        in >> rocks[i].mass >> rocks[i].energy; //mass and energy for each rock
        //rocks[i].printRock(&out);
    }

    int rocksPerCity;
    int city;
    cities = new City[nCitiesTot];

    for (int i = 0; i < nDifferentRocks; ++i) { //i-th rock
        in >> rocksPerCity;

        //out << rocksPerCity << endl;
        for (int j = 0; j < rocksPerCity; ++j) {
            in >> city;

            cities[city].rocks.push_back(rocks[i]); //push(number of rock)
            //out << city << " ";
        }
    }

    weights = new int *[nCitiesTot - 1]; //1 is 0

    for (int i = 1; i < nCitiesTot; ++i) {
        weights[i] = new int[i];
        for (int j = 0; j < i; ++j) {
            in >> weights[i][j];
            //out << weights[i][j] << " ";
        }
        //out << endl;
    }
}

//@todo adapt to new structures
void output(double E, double G, double T, int *rocks, int *cities) {
    ofstream out("output.txt");

    //final glove's energy, rocks' energy and time for the journey
    out << scientific << setprecision(10) << E << " ";
    out << scientific << setprecision(10) << G << " ";
    out << scientific << setprecision(10) << T << endl;

    //print all taken rocks
    for (int i = 0; i < nTakenRocks; i++) {
        printf("%d ", rocks[i]);
    }
    printf("\n");

    //print the path (ordered, t0 = tN = source)
    for (int i = 0; i < nCitiesTot; i++) {
        printf("%d ", cities[i]);
    }
    printf("\n***\n");
}

int calcV(int vMax, int vMin, int W, int C) {
    return vmax - W * ((vMax - vMin) / C);
}

int main() {
    int *rocks, cities;
    input();


    return 0;
}

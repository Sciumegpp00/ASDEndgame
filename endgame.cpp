#include <iostream>
#include <stdlib.h>
#include "endgame.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>

#define INPUT_FILE "../inputs/input6.txt"
#define OUTPUT_FILE "../outputs/output6.txt"
//#define INPUT_FILE "input.txt"
//#define OUTPUT_FILE "output.txt"

using namespace std;

struct Rock {
    int mass;
    int energy;
    int availability;
    float ratio;
    short id;

    //variable for rapport velocity/distance (depends from path)

    /*Rock(int m, int e){
        mass = m;
        energy = e;
        k = energy/mass;
    }*/

    void printRock(ofstream *out) {
        *out << mass << " " << energy << endl;
    }
};

struct City {
    vector<Rock> rocks;
};

struct CitySol {
    short city;
    short rock;
};

//struct Edge {
//    int weight;
//    short c1;
//    short c2;
//};

int getDistance(short city1, short city2);
void input();
void output(CitySol* solution, short* orderedRocks, double gloveEnergy, double rocksEnergy, double time);
void bbTsp(vector<short>* path, int cost, vector<short>* remaining, int n, int i);
void removeVector(vector<short>* v, short x);
int simpleMinPath(CitySol* sol, short* rocksSol);
void printMinPath(vector<short>* minSol);
void printMinCost(int minCost);
void printSolution(CitySol* solution, short* orderedRocks, double gloveEnergy, double rocksEnergy, double time);
double calcSpeed(int W);
double calcEnergy(int rocksEnergy, double totalTime);
void findRocks(CitySol* city, vector<Rock>* cityRocks, short* rocksSol, int* capacityLeft);
void calcFinalResources(CitySol* sol, double* totalTime, double* totalEnergy, int* rocksEnergy);


short nCitiesTot;
short nDifferentRocks;
short source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
//TODO: idea: CitySol *solution che contiene in ordine di id citta (seguendo percorso minimo) l'id città e la roccia presa (o non presa) --> quindi è gran parte della soluzione finale
int **weights;
//int globalMinWeight;
//int minCost; //costo di un percorso minimo a caso (il primo percorso che prendo)
//vector<short>* minSol = nullptr;
//vector<short>* remainingBbTsp = nullptr;


int main() {
    CitySol* solution;
    short* rocksSolution;
    int cost, rocksEnergy;
    double totalTime, totalEnergy;

    input();

    solution = new CitySol[nCitiesTot+1]; // Last city is the source
    rocksSolution = new short[nDifferentRocks];

    cost = simpleMinPath(solution, rocksSolution);
    calcFinalResources(solution, &totalTime, &totalEnergy, &rocksEnergy);
    //bbTsp(&sol, 0, remainingBbTsp, nCitiesTot-1, 1);

    cout << "Costo: " << cost << endl;
    printSolution(solution, rocksSolution, totalEnergy, rocksEnergy, totalTime);
    output(solution, rocksSolution, totalEnergy, rocksEnergy, totalTime);

    //bestRocksForPath(minSol);

    //output();
    return 0;
}

void input() {
    ifstream in(INPUT_FILE);

    in >> nCitiesTot >> source;
    in >> nDifferentRocks >> C >> R >> vmin >> vmax;

    rocks = new Rock[nDifferentRocks];

    for (int i = 0; i < nDifferentRocks; ++i) {
        rocks[i].id = i;
        in >> rocks[i].mass >> rocks[i].energy; //mass and energy for each rock
        rocks[i].ratio = (float) rocks[i].energy/rocks[i].mass;
//        cout << "rock n " << rocks[i].id << " mass: " << rocks[i].mass << " energy: " << rocks[i].energy << " ratio: " << rocks[i].ratio << endl;
    }

    int city;
    cities = new City[nCitiesTot];

    for (int i = 0; i < nDifferentRocks; ++i) { //i-th rock
        in >> rocks[i].availability;

        for (int j = 0; j < rocks[i].availability; ++j) {
            in >> city;
            cities[city].rocks.push_back(rocks[i]);
        }
    }

    weights = new int *[nCitiesTot - 1]; //1 is 0
//    remainingBbTsp = new vector<short>();

    //globalMinWeight = INT_MAX;

    //inserimento dei pesi tra città
    for(short i=0; i<nCitiesTot-1; i++){
        weights[i] = new int[i+1];
        for(short j=0; j<i+1; j++){
            in >> weights[i][j];

//            if(weights[i][j] < globalMinWeight)
//                globalMinWeight = weights[i][j];
        }
    }
}

void output(CitySol* solution, short* orderedRocks, double gloveEnergy, double rocksEnergy, double time){
    ofstream out(OUTPUT_FILE);

//    //final glove's energy, rocks' energy and time for the journey
    out << scientific << setprecision(10) << gloveEnergy << " "; //total glove's energy
    out << scientific << setprecision(10) << rocksEnergy << " "; //total rocks energy
    out << scientific << setprecision(10) << time << endl; //time used

    for (int i = 0; i < nDifferentRocks; ++i) {
        out << orderedRocks[i] << " ";
    }

    out << endl;

    for (int i = 0; i <= nCitiesTot; ++i) {
        out << solution[i].city << " ";
    }

    out << "\n***";
}

void printSolution(CitySol* solution, short* orderedRocks, double gloveEnergy, double rocksEnergy, double time){
//    ofstream out(OUTPUT_FILE);

//    //final glove's energy, rocks' energy and time for the journey
//    cout << scientific << setprecision(10) << gloveEnergy << " "; //total glove's energy
//    cout << scientific << setprecision(10) << rocksEnergy << " "; //total rocks energy
//    cout << scientific << setprecision(10) << time << endl; //time used

    cout << "Total energy: " << gloveEnergy << endl;
    cout << "Rocks energy: " << rocksEnergy << endl;
    cout << "Time: " << time << endl;

    cout << "Rocce: " << nDifferentRocks << endl;
    for (int i = 0; i < nDifferentRocks; ++i) {
        cout << orderedRocks[i] << " ";
    }

    cout << "\nPercorso: " << nCitiesTot << endl;
    for (int i = 0; i <= nCitiesTot; ++i) {
        cout << solution[i].city << " ";
    }

    cout << "\n***";
}

int simpleMinPath(CitySol* sol, short* rocksSol) { //int rocksWeight, int C
    short i, j, swap, minCityIndex;
    int minDistance, localDistance, cost = 0, capacityLeft = C;
//    double speed, energy, time;

    for(i = 0; i < nCitiesTot; i++) {
        sol[i].city = i;
    }

    sol[0].city = source;
    sol[source].city = 0;
    sol[nCitiesTot].city = source;

    for (i = 0; i < nDifferentRocks; ++i) {
        rocksSol[i] = -1;
    }

    for(i = nCitiesTot-1; i >= 0; i--) {
        minCityIndex = i;
        minDistance = getDistance(sol[i + 1].city, sol[i].city); //minTime = getDistance(sol[i+1].city, sol[i].city)/calcSpeed(W);

        for(j = i-1; j > 0; j--) {
            localDistance = getDistance(sol[i + 1].city, sol[j].city);
            if(localDistance < minDistance) {
                minCityIndex = j;
                minDistance = localDistance;
            }
        }
        cost += minDistance;

        swap = sol[i].city;
        sol[i].city = sol[minCityIndex].city;
        sol[minCityIndex].city = swap;

        // Take rocks
        findRocks(&sol[i], &(cities[sol[i].city].rocks), rocksSol, &capacityLeft);
    }

    return cost;
}

void findRocks(CitySol* city, vector<Rock>* cityRocks, short* rocksSol, int* capacityLeft) {
    if(cityRocks->empty()) {
        city->rock = -1;
        return;
    }

    cityRocks = &(cities[city->city].rocks);

    // Find maximum ratio for all rock of the city
    Rock maxRock = {};
    maxRock.ratio = 0;

    for (auto & cityRock : *cityRocks) {
        if(*capacityLeft - cityRock.mass < 0 || rocksSol[cityRock.id] != -1)
            continue;
        if(cityRock.ratio > maxRock.ratio)
            maxRock = cityRock;
    }

    if(maxRock.ratio == 0)
        city->rock = -1;
    else {
        city->rock = maxRock.id;
        rocksSol[city->rock] = city->city;
        *capacityLeft -= rocks[city->rock].mass;
    }
}

void calcFinalResources(CitySol* sol, double* totalTime, double* totalEnergy, int* rocksEnergy) {
    int rocksWeight = 0, distance;
    double speed, time;

    *rocksEnergy = 0;

    for (short i = 0; i < nCitiesTot; i++) {
        if(sol[i].rock != -1) {
            rocksWeight += rocks[sol[i].rock].mass;
            *rocksEnergy += rocks[sol[i].rock].energy;
        }

        distance = getDistance(sol[i].city, sol[i+1].city);
        speed = calcSpeed(rocksWeight);
        time = distance / speed;

        *totalTime += time;
    }

    *totalEnergy = calcEnergy(*rocksEnergy, *totalTime);
}


void simpleMinPath(vector<CitySol>* sol, vector<short>* remaining) { //first call of function: sol with only source inside and remaining with all but source

    //minimum path for each city
    /*if(sol->size() < nCitiesTot){
        CitySol minCity;
        short presentCity = sol->back().city;
        int presentWeight, minWeight;
        if(sol->back().city == sol->front().city && sol->back().city == 0){
            minWeight = INT16_MAX;
        }else{
            minWeight = getWeight(sol->back().city, remaining->back());
            minCity.city = remaining->back();
        }

        for (auto it = remaining->begin() ; it != remaining->end(); ++it) {
            presentWeight = getDistance(*it, presentCity);
            if(presentWeight < minWeight){
                minCity.city = *it;
                minWeight = presentWeight;
            }
        }*/

//        sol->push_back(minCity);
//        removeVector(remaining, minCity.city);
//        minCost += minWeight;
//        simpleMinPath(sol, remaining);
//    }

}


int getDistance(short city1, short city2) {
    return city1 < city2
           ? weights[city2 - 1][city1]
           : weights[city1 - 1][city2];
}

double calcSpeed(int W) {
    if(C != 0)
        return vmax - W * ((vmax - vmin) / C);
    return vmax;
}

double calcEnergy(int rocksEnergy, double totalTime) {
    return rocksEnergy - R * totalTime;
}

void removeVector(vector<short>* v, short x) {
    for (auto it = v->begin() ; it != v->end(); it++) {
        if(*it == x) {
            v->erase(it);
            return;
        }
    }
}

void printMinPath(vector<short>* minSol){
    cout << "\nPercorso minimo: " << endl;
    for(int i = 0; i < minSol->size(); i++){
        cout << (*minSol)[i] << " ";
    }
    cout << source << endl;
}





int calcLb(short origin, vector<short>* remaining, int cost) {
    int out, back, transfer, costLocal, minTransferOut, localTransferOut;

    //out = INT_MAX;
    //back = INT_MAX;
//    transfer = INT_MAX;
//    out = getDistance(origin, (*remaining)[0]);
//    back = getDistance((*remaining)[0], source);
//    transfer = out + getDistance((*remainingBbTsp)[0], (*remainingBbTsp)[1]);

    for (int i = 0; i < remaining->size(); i++) {
        // calculate minimum cost of exiting edges at this point
        costLocal = getDistance(origin, (*remaining)[i]);
        if(costLocal < out) {
            out = costLocal;
        }

        // Calculate for each city the best way to pass into it (best duo entry + exit)
        // origin --> qualsiasi nodo + qualsiasi nodo
//        minTransferOut = INT_MAX;
//        for (int j = 0; j < remaining->size(); j++){ //consider the remainingBbTsp nodes
//            if(i != j){
//                localTransferOut = getDistance((*remaining)[i], (*remaining)[j]);
//                if(localTransferOut < minTransferOut)
//                    minTransferOut = localTransferOut;
//            }
//        }
//        transfer += minTransferOut;

        //calculate lb cost to go back to source of the nearest cities
        costLocal = getDistance((*remaining)[i], source);
        if(costLocal < back) {
            back = costLocal;
        }
    }

//    transfer = (globalMinWeight * remaining->size());
    // TODO: Try to remove /2 for remainingBbTsp nodes > 10 or 15 for example
//    return cost + ((out + back + transfer) / 2);
    return cost + out + back + transfer;
}

void bbTsp(vector<short>* path, int cost, vector<short>* remaining, int n, int i) {
//    short origin = path->back();
//    vector<short> choices(*remaining);
//    int lb, localCost;
//
//    for (auto it = choices.begin() ; it != choices.end(); ++it) {
//        localCost = cost + getDistance(origin, *it);
//        path->push_back(*it);
//        removeVector(remaining, *it);
//
//        if(i < n) {
//            lb = calcLb(*it, remaining, localCost);
//            if(lb < minCost){
//                bbTsp(path, localCost, remaining, n, i + 1);
//            }
//        } else {
//            cost = localCost;
//            cost += getDistance(*it, (*path)[0]); //path->front() Returns a reference to the first element in the vector.
//            if(cost < minCost) {
//                if(minSol) {
//                    delete minSol;
//                    minSol = nullptr;
//                }
//                minSol = new vector<short>(*path);
//                minCost = cost;
//            }
//        }
//
//        path->pop_back();
//        remaining->push_back(*it);
//    }
}

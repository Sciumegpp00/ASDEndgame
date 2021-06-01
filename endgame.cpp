#include <iostream>
#include <cstdlib>
#include "endgame.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <queue>

#ifdef DEBUG
#define INPUT_FILE "../inputs/input9.txt"
#define OUTPUT_FILE "../outputs/output9.txt"
#endif // DEBUG

#ifndef DEBUG
#define INPUT_FILE "input.txt"
#define OUTPUT_FILE "output.txt"
#endif // DEBUG

using namespace std;

struct Rock {
    int mass;
    int energy;
    float ratio;
    bool toTake = false;
    short availability;
    short *cities;
};

struct RockRatio {
    short id;
    float ratio;
};

bool operator>(const Rock &r1, const Rock &r2);

struct City {
    short nRocks = 0;
    short *rocks;
    short addedRocks = 0;
};

struct CitySol {
    short cityId;
    short rockId;
};


int getDistance(short city1, short city2);
void input();
void output(CitySol *solution, short *orderedRocks, double gloveEnergy, double rocksEnergy, double time);
void bbTsp(vector<short> *path, int cost, vector<short> *remaining, int n, int i);
void removeVector(vector<short> *v, short x);
int simpleMinPath(CitySol *sol, short *rocksSol);
void printMinPath(vector<short> *minSol);
void printMinCost(int minCost);
void printSolution(CitySol *solution, short *orderedRocks, double gloveEnergy, double rocksEnergy, double time);
double calcSpeed(int W);
double calcEnergy(int rocksEnergy, double totalTime);
void findRocks(CitySol *citySol, short *rocksSol, int *capacityLeft);
void calcFinalResources(CitySol *sol, double *totalTime, double *totalEnergy, int *rocksEnergy);
void chooseRocksToTake();


short nCitiesTot;
short nDifferentRocks;
short source;
int C;
double R;
double vmax;
double vmin;
Rock *rocks;
City *cities;
int **weights;
RockRatio *rocksRatio;
//int globalMinWeight;
//int minCost; //costo di un percorso minimo a caso (il primo percorso che prendo)
//vector<short>* minSol = nullptr;
//vector<short>* remainingBbTsp = nullptr;


//bool operator>(const Rock &r1, const Rock &r2) {
//    return r1.ratio < r2.ratio;
//}

int minRockRatioComparator(const RockRatio *r1, const RockRatio *r2) {
    return r1->ratio < r2->ratio
           ? -1
           : 1;
}



int main() {
//    Rock test[5];
//
//    test[0].id = 0; test[0].ratio = 4.6;
//    test[1].id = 1; test[1].ratio = 2.3;
//    test[2].id = 2; test[2].ratio = 10.55;
//    test[3].id = 3; test[3].ratio = 4.9;
//    test[4].id = 4; test[4].ratio = 6.8;
//
//    for (int i = 0; i < 5; ++i) {
//        cout << test[i].ratio << " ";
//    }
//    cout << endl;
//
//    qsort(test, 10, sizeof(Rock),
//          reinterpret_cast<int (*)(const void *, const void *)>(minRockRatioComparator));

    CitySol *solution;
    short *rocksSolution;
    int cost, rocksEnergy;
    double totalTime = 0, totalEnergy = 0;
#ifdef DEBUG
    clock_t startTime, afterInput, afterSolution, afterFinalCalc;
    startTime = clock();
#endif // DEBUG

    input();
    chooseRocksToTake();

#ifdef DEBUG
    afterInput = clock();
#endif // DEBUG

    solution = new CitySol[nCitiesTot + 1]; // Last cityId is the source
    rocksSolution = new short[nDifferentRocks];
    cost = simpleMinPath(solution, rocksSolution);

#ifdef DEBUG
    afterSolution = clock();
#endif // DEBUG

    calcFinalResources(solution, &totalTime, &totalEnergy, &rocksEnergy);
    //bbTsp(&sol, 0, remainingBbTsp, nCitiesTot-1, 1);

#ifdef DEBUG
    afterFinalCalc = clock();
#endif // DEBUG

    output(solution, rocksSolution, totalEnergy, rocksEnergy, totalTime);

#ifdef DEBUG
    cout << "Costo: " << cost << endl;
    printSolution(solution, rocksSolution, totalEnergy, rocksEnergy, totalTime);
    cout << endl << endl;
    cout << "Input duration: " << afterInput - startTime << endl;
    cout << "Solution creation duration: " << afterSolution - afterInput << endl;
    cout << "Solution calc duration: " << afterFinalCalc - afterSolution << endl;
#endif // DEBUG

    return 0;

}

void input() {
#ifdef DEBUG
    clock_t startTime = clock(), afterRocksInput, afterRocksPerCity, afterRocksSort, afterWeightInput;
#endif // DEBUG

    ifstream in(INPUT_FILE);

    in >> nCitiesTot >> source;
    in >> nDifferentRocks >> C >> R >> vmin >> vmax;

    rocks = new Rock[nDifferentRocks];
    rocksRatio = new RockRatio[nDifferentRocks];

    for (short i = 0; i < nDifferentRocks; ++i) {
        in >> rocks[i].mass >> rocks[i].energy; //mass and energy for each rockId
        rocks[i].ratio = (float) rocks[i].mass / (float) rocks[i].energy;
        //vettore rocce id ratio
        rocksRatio[i].id = i;
        rocksRatio[i].ratio = rocks[i].ratio;
//        cout << "rockId n " << rocks[i].id << " mass: " << rocks[i].mass << " energy: " << rocks[i].energy << " ratio: " << rocks[i].ratio << endl;
    }

#ifdef DEBUG
    for (int i = 0; i < nDifferentRocks; ++i) {
        cout << rocksRatio[i].ratio << " ";
    }

    cout << endl;
    afterRocksInput = clock();
#endif // DEBUG

    short city;
    cities = new City[nCitiesTot];

//    for (int i = 0; i < nDifferentRocks; ++i) { //i-th rockId
//        in >> rocks[i].availability;
//
//        for (int j = 0; j < rocks[i].availability; ++j) {
//            in >> city;
//        }
//    }

    for (short i = 0; i < nDifferentRocks; ++i) { //i-th rockId
        in >> rocks[i].availability;
        rocks[i].cities = new short[rocks[i].availability];

        for (short j = 0; j < rocks[i].availability; ++j) {
            in >> city;
            rocks[i].cities[j] = city;
            cities[city].nRocks++;
        }
    }

    for (short i = 0; i < nCitiesTot; ++i) {
        cities[i].rocks = new short[cities[i].nRocks];
    }

    for (short i = 0; i < nDifferentRocks; ++i) {
        for (short j = 0; j < rocks[i].availability; ++j) {
            city = rocks[i].cities[j];
            cities[city].rocks[cities[city].addedRocks++] = i;
        }
        delete[] rocks[i].cities;
    }

#ifdef DEBUG
    afterRocksPerCity = clock();
#endif // DEBUG

//    for (int i = 0; i < nCitiesTot; ++i) {
//        qsort(cities[i].rocks, cities[i].nRocks, sizeof(short),
//              reinterpret_cast<int (*)(const void *, const void *)>(minRockRatioComparator));
//    }

#ifdef DEBUG
    afterRocksSort = clock();
#endif // DEBUG

    weights = new int *[nCitiesTot - 1]; //1 is 0

    //inserimento dei pesi tra città
    for (short i = 0; i < nCitiesTot - 1; i++) {
        weights[i] = new int[i + 1];
        for (short j = 0; j < i + 1; j++) {
            in >> weights[i][j];

//            if(weights[i][j] < globalMinWeight)
//                globalMinWeight = weights[i][j];
        }
    }

#ifdef DEBUG
    afterWeightInput = clock();

    cout << "Input rocks: " << afterRocksInput - startTime << endl;
    cout << "Rocks per cityId: " << afterRocksPerCity - afterRocksInput << endl;
    cout << "Rock sort: " << afterRocksSort - afterRocksPerCity << endl;
    cout << "Weight input: " << afterWeightInput - afterRocksSort << endl;
    cout << endl << endl;
#endif // DEBUG
}

void output(CitySol *solution, short *orderedRocks, double gloveEnergy, double rocksEnergy, double time) {
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
        out << solution[i].cityId << " ";
    }

    out << "\n***";
}

void printSolution(CitySol *solution, short *orderedRocks, double gloveEnergy, double rocksEnergy, double time) {
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
        cout << solution[i].cityId << " ";
    }

    cout << "\n***";
}

void chooseRocksToTake() {
    int capacityLeft = C, taken = 0;

    qsort(rocksRatio, nDifferentRocks, sizeof(RockRatio),
          reinterpret_cast<int (*)(const void *, const void *)>(minRockRatioComparator));

    for (int i = 0; i < nDifferentRocks && capacityLeft != 0; ++i) {
        if(rocks[rocksRatio[i].id].mass > capacityLeft)
            continue;

        capacityLeft -= rocks[rocksRatio[i].id].mass;
        rocks[rocksRatio[i].id].toTake = true;
        taken++;
    }

    cout << "Taken rocks: " << taken << endl;
}

int simpleMinPath(CitySol *sol, short *rocksSol) { //int rocksWeight, int C
    short i, j, swap, minCityIndex;
    int minDistance, localDistance, cost = 0, capacityLeft = C;
//    double speed, energy, time;

    for (i = 0; i < nCitiesTot; i++) {
        sol[i].cityId = i;
    }

    sol[0].cityId = source;
    sol[source].cityId = 0;
    sol[nCitiesTot].cityId = source;

    for (i = 0; i < nDifferentRocks; ++i) {
        rocksSol[i] = -1;
    }

    for (i = nCitiesTot - 1; i >= 0; i--) {
        minCityIndex = i;
        minDistance = getDistance(sol[i + 1].cityId, sol[i].cityId); //minTime = getDistance(sol[i+1].cityId, sol[i].cityId)/calcSpeed(W);

        for (j = i - 1; j > 0; j--) {
            localDistance = getDistance(sol[i + 1].cityId, sol[j].cityId);
            if (localDistance < minDistance) {
                minCityIndex = j;
                minDistance = localDistance;
            }
        }
        cost += minDistance;

        swap = sol[i].cityId;
        sol[i].cityId = sol[minCityIndex].cityId;
        sol[minCityIndex].cityId = swap;

        // Take rocks
        findRocks(&sol[i], rocksSol, &capacityLeft);
    }

    return cost;
}

void findRocks(CitySol *citySol, short *rocksSol, int *capacityLeft) { //vector<Rock>* cityRocks
    auto cityRocks = cities[citySol->cityId].rocks;
    auto nRocks = cities[citySol->cityId].nRocks;

    if (nRocks == 0 || *capacityLeft == 0) { // cityRocks->empty()
        citySol->rockId = -1;
        return;
    }

    // Find maximum ratio for all rockId of the cityId

//    for (int i = 0; i < nRocks; ++i) {
//        if ((*capacityLeft - rocks[cityRocks[i]].mass) < 0 || rocksSol[cityRocks[i]] != -1)
//            continue;
//
//        citySol->rockId = cityRocks[i];
//        rocksSol[citySol->rockId] = citySol->cityId;
//        *capacityLeft -= rocks[citySol->rockId].mass;
//        return;
//    }
//    citySol->rockId = -1;

//    for (int i = 0; i < cityRocks->size(); ++i) {
//        if (*capacityLeft - ((*cityRocks).top()).mass < 0 || rocksSol[(*cityRocks).top().id] != -1)
//            continue;
//
//        cityId->rockId = (*cityRocks).top().id;
//        rocksSol[cityId->rockId] = cityId->cityId;
//        *capacityLeft -= rocks[cityId->rockId].mass;
//        return;
//    }


    Rock* localRock;
    float localRatio;
    short maxRockId;
    float maxRockRatio = 0;

    for (int i = 0; i < nRocks; ++i){
        localRock = rocks + cityRocks[i];
//        if(*capacityLeft - rocks[cityRocks[i]].mass < 0 || rocksSol[cityRocks[i]] != -1)
        if(!localRock->toTake)
            continue;

        if(localRock->availability < 20) {
            maxRockId = cityRocks[i];
            maxRockRatio = localRatio;
            break;
        }

        localRatio = localRock->ratio / (float) localRock->availability;
        localRock->availability--;
        if(localRatio > maxRockRatio) {
            maxRockId = cityRocks[i];
            maxRockRatio = localRatio;
        }
    }

    if(maxRockRatio == 0)
        citySol->rockId = -1;
    else {
        citySol->rockId = maxRockId;
        rocksSol[maxRockId] = citySol->cityId;
        *capacityLeft -= rocks[maxRockId].mass;
        rocks[maxRockId].toTake = false;
    }
}

//int fullMinPathSearch(CitySol* sol, short* rocksSol) {
//
//}

void calcFinalResources(CitySol *sol, double *totalTime, double *totalEnergy, int *rocksEnergy) {
    int rocksWeight = 0, distance;
    double speed, time;

    *rocksEnergy = 0;
    *totalTime = 0;

    for (short i = 0; i < nCitiesTot; i++) {
        if (sol[i].rockId != -1) {
            rocksWeight += rocks[sol[i].rockId].mass;
            *rocksEnergy += rocks[sol[i].rockId].energy;
        }

        distance = getDistance(sol[i].cityId, sol[i + 1].cityId);
        speed = calcSpeed(rocksWeight);
        time = distance / speed;

        *totalTime += time;
    }

    *totalEnergy = calcEnergy(*rocksEnergy, *totalTime);
}


int getDistance(short city1, short city2) {
    return city1 < city2
           ? weights[city2 - 1][city1]
           : weights[city1 - 1][city2];
}

double calcSpeed(int W) {
    if (C != 0)
        return vmax - W * ((vmax - vmin) / C);
    return vmax;
}

double calcEnergy(int rocksEnergy, double totalTime) {
    return rocksEnergy - R * totalTime;
}

void removeVector(vector<short> *v, short x) {
    for (auto it = v->begin(); it != v->end(); it++) {
        if (*it == x) {
            v->erase(it);
            return;
        }
    }
}


int calcLb(short origin, vector<short> *remaining, int cost) {
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
        if (costLocal < out) {
            out = costLocal;
        }

        // Calculate for each cityId the best way to pass into it (best duo entry + exit)
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
        if (costLocal < back) {
            back = costLocal;
        }
    }

//    transfer = (globalMinWeight * remaining->size());
    // TODO: Try to remove /2 for remainingBbTsp nodes > 10 or 15 for example
//    return cost + ((out + back + transfer) / 2);
    return cost + out + back + transfer;
}

void bbTsp(vector<short> *path, int cost, vector<short> *remaining, int n, int i) {
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

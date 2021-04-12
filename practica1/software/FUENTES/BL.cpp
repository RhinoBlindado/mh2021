// INTRODUCTION

/**
 *  [CASTELLANO]
 * 
 * Practica 1.b - Algoritmos Greedy y Busqueda Local para PAR
 * Asignatura: Metaheuristicas
 * Autor: Valentino Lugli (Github: @RhinoBlindado)
 * Marzo, Abril 2021
 */

/**
 *  [ENGLISH]
 * 
 * Practice 1.b - Greedy and Local Search Algorithms for Clustering with Restrictions
 * Course: Metaheuristics
 * Author: Valentino Lugli (Github: @RhinoBlindado)
 * March, April 2021
 */

// LIBRARIES
//  I/O.
#include <iostream>
//  Vectors.
#include <vector>   
//  Open files.
#include <fstream>
//  Use StringSteam to parse data.  
#include <sstream>
//  Randomness.
#include "libs/random.h"
//  Take the time between executions.
#include <ctime>
//  Do the random shuffle.
#include <algorithm> 
//  Get max size of data types.
#include <limits>
//  Perform math calculations.
#include <cmath>
//  Use the Pair class
#include <utility>

using namespace std;

//  STRUCTS
/**
 * @brief Small container for the restrictions
 * @param   x_0     ith Row in the Restriction Matrix
 * @param   x_1     jth Row in the Restriction Matrix
 * @param   R       The kind of restriction held 
 */
struct triplet
{
    int x_0;
    int x_1;
    int R;
};

// FUNCTIONS
/**
 * @brief Parse and fill the matrices with the supplied data.
 * @param dataSetPath               Path that leads to a *.set file.
 * @param dataSetRestrictiosPath    Path that leads to a *.const file.
 * @param X                         Matrix that will be filled with n instances of d dimensions.
 * @param MR                        Matrxi that will be filled with n * n instances with their restrictions.
 */
void loadData(string dataSetPath, string dataSetRestrictionsPath, vector<vector<float>> &X, vector<vector<int>> &MR)
{
    ifstream dataSetFile (dataSetPath);
    ifstream dataSetRestrictions (dataSetRestrictionsPath);
    string line, item;

    if(dataSetFile.is_open())
    {
        int i=0, j=0;
        while(getline(dataSetFile, line, '\n'))
        {
            stringstream aux(line);
            while (getline(aux, item, ','))
            {
                X[i][j] = stof(item);
                j++;
            }
            i++;
            j=0;
        }
    }
    else
    {
        cout<<"Error: "<<dataSetPath<<" is not a valid file."<<endl;
        exit(-1);
    }

    if(dataSetRestrictions.is_open())
    {
        int i=0, j=0;
        while(getline(dataSetRestrictions, line, '\n'))
        {
            stringstream aux(line);
            while (getline(aux, item, ','))
            {
                MR[i][j] = stof(item);
                j++;
            }
            i++;
            j=0;
        }
    }
    else
    {
        cout<<"Error: "<<dataSetRestrictionsPath<<" is not a valid file."<<endl;
        exit(-1);
    }
}

/**
 * @brief Fill a list containing all the restrictions in the form of (x_i, x_j, R), also one with only MUST LINK and CANNOT LINK restrictions.
 * @param n     Number of instances.
 * @param MR    Restriction matrix.
 * @param LR    Vector of Triplets to be filled.
 * @param ML    Vector of MUST LINK Triplets to be filled. 
 * @param CL    Vector of CANNOT LINK Triplets to be filled.
 */
void fillTriplet(int n, vector<vector<int>> MR, vector<triplet> &LR, vector<triplet> &ML, vector<triplet> &CL)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            if(i != j && MR[i][j] != 0)
            {
                triplet aux = {i, j, MR[i][j]};
                LR.push_back(aux);
                if(MR[i][j] == 1)
                {
                    ML.push_back(aux);
                }

                if(MR[i][j] == -1)
                {
                    CL.push_back(aux);
                }
            }
        }
    }
}

/**
 * @brief Generate the lambda component of the fitness calculation.
 * @param n     Number of instances.
 * @param dim   Number of dimenions.
 * @param X     List of instances with their d dimensions.
 * @param rSize Total number of restrictions.
 * @return Lambda, calculated as D, the largest distance in the current dataset divided by the number of restrictions.
 */
double genLambda(int n, int dim, vector<vector<float>> X, int rSize)
{
    double maxDist = -1;
    double actDist;

    /*
     * Lambda = D / |R|
     *  - D is the largest distance between two points.
     *  - |R| is the number of restrictions.
     */

    // Calculating maxDist 
    for (int i = 0; i < n; i++)
    {
        for(int j = i; j < n; j++)
        {

            if(i != j)
            {
                actDist = 0;
                for(int k = 0; k < dim; k++)
                {
                    actDist += pow(X[i][k] - X[j][k], 2);
                }

                actDist = sqrt(actDist);

                if(actDist > maxDist)
                {
                    maxDist = actDist;
                }
            }
        }
    }

    maxDist = ceil(maxDist);

    return maxDist / (double) rSize;
    
}

/**
 * @brief Random Wrapper for the RNG Library to be used in the Random Shuffle function
 * @param i     Index that is used for Random Shuffle to generate a number. 
 * @return Random integer between 0 and i
 */
int rWrapper(int i)
{
    return Randint(0, i);
}

/**
 * @brief Generate the Initial Solution
 * @param S             Vector with size n (wih n >= numClusters)
 * @param numClusters   Number of Clusters (k)
 * @param centroidList  List that maitains what instances are assigned to a cluster
 */
void genInitSolution(vector<int> &S, int numClusters, vector<int> &centroidList)
{
    int size = S.size();

    // Fill the S vector with the ith cluster so that at least each of the clusters appears once
    for (int i = 0; i < numClusters; i++)
    {
        S[i] = i;
        centroidList[i]++;
    }

    // Fill the rest of the vector with a random cluster
    for(int i = numClusters; i < size; i++)
    {
        S[i] = Randint(0, numClusters-1);
        centroidList[S[i]]++;
    }
    
    // Shuffle the vector to randomize it
    random_shuffle(S.begin(), S.end(), rWrapper);
   
}

/**
 * @brief Generate a "virtual" neigborhood of the current solution
 * @param S             Vector of the current solution.
 * @param sizeK         Number of centroids.
 * @param clusterCount  Vector that mantains the number of instances assigned to each cluster.
 * @return  A vector of pairs [l, i]; a change to S in the form of S[l] = i.
 */
vector<pair<int,int>> genNeighbor(vector<int> S, int sizeK, vector<int> clusterCount)
{
    vector<pair<int,int>> neigh;
    pair<int,int> aux;

    for (int i = 0; i < S.size(); i++)
    {
        for (int j = 0; j < sizeK; j++)
        {
            /** 
             * For each position of S, generate a neighbor, that is: changing that instance from a cluster to every other one.
             * Except if that cluster only has 1 instance; then do nothing because all clusters must have at least one instance assigned.
             */
            if(S[i] != j && clusterCount[S[i]] > 1)
            {
                aux.first = i;
                aux.second = j;
                neigh.push_back(aux);
            }
            
        }
        
    }

    // Shuffle the resulting vector.
    random_shuffle(neigh.begin(), neigh.end(), rWrapper);
    
    return neigh;
}

/**
 * @brief Calculates the infeasibility of the current solution.
 * @param S     Vector of the current solution.
 * @param ML    Vector of MUST-LINK Triplets.
 * @param CL    Vector of CANNOT-LINK Triplets.
 * @return  The count of how many times a restrictions has been violated.
 */
int infeasibility(vector<int> S, vector<triplet> ML, vector<triplet> CL)
{
    int count = 0;
    int mlSize = ML.size(),
        clSize = CL.size();

    // Checking MUST-LINK restrictions.
    for(int i = 0; i < mlSize; i++)
    {
        // If x_0 and x_1 belong to different clusters, ML is violated.
        if(S[ML[i].x_0] != S[ML[i].x_1])
        {
            count++;
        }
    }

    // Checking CANNOT-LINK restrictions.
    for(int i = 0; i < clSize; i++)
    {
        // If x_0 and x_1 belong to the same cluster, CL is violated.
        if(S[CL[i].x_0] == S[CL[i].x_1])
        {
            count++;
        }

    }

    return count;
}


/**
 * @brief Calculate the instances associated with a centroid given the current solution.
 * @param actInst   The vecto f actual instances associated with the centroid to be filled.
 * @param X         Array that contains n instances with their d dimensions. 
 * @param actK      The current centroid/cluster
 * @param S         Vector that contains the actual solution
 */
void calcInstances(vector<int> &actInst, vector<vector<float>> X, int actK, vector<int> S)
{
    vector<int> actList;
    int size = S.size();

    for (int i = 0; i < size; i++)
    {
        if(S[i] == actK)
        {
            actInst.push_back(i);
        }
    }

}

/**
 * @brief Calculate the coordinates of a centroid given the list of instances associated with it.
 * @param actCoord  The centroid coordinates to be calculated.
 * @param actInst   The actual instances associated with the centroid.
 * @param X         Array that contains n instances with their d dimensions. 
 */
void calcCentroidCoords(vector<float> &actCoord, vector<int> &actInst, vector<vector<float>> X)
{
    int dimensions = X[0].size();
    int instances = actInst.size();

    float sum;
    
    for (int i = 0; i < dimensions; i++)
    {
        sum = 0;

        for(int j = 0; j < instances; j++)
        {
            sum += X[actInst[j]][i];
        }

        actCoord[i] = sum / instances;
    }

}

/**
 * @brief Get the average intracluster distance of a certain cluster.
 * @param actCoord  The actual centroid coordinates.
 * @param actInst   The actual instances associated with the centroid.
 * @param X         Array that contains n instances with their d dimensions. 
 * @return Average intracluster distane of a certain cluster.
 */
double calcIntraDiff(vector<float> actCoord, vector<int> actInst, vector<vector<float>> X)
{
    double sum = 0;
    double actDist;

    int dimensions = X[0].size();
    int instances = actInst.size();

    for(int i = 0; i < instances; i++)
    {
        actDist = 0;
        for(int j = 0; j < dimensions; j++)
        {
            actDist += pow(X[actInst[i]][j] - actCoord[j], 2);
        }
        actDist = sqrt(actDist);
        sum += actDist;
    }

    return sum / instances;
}

/**
 * @brief Get the average intracluster distance of the current solution.
 * @param S         Vector that contains the actual solution
 * @param X         Array that contains n instances with their d dimensions. 
 * @param k         Number of Centroids.
 * @return Average intracluster distance of the current solution.
 */
double intraClusterDistance(vector<int> S, vector<vector<float>> X, int k)
{
    vector<int> actInstance;
    vector<float> actCoords (X[0].size(), 0);
    double  actDiff,
            totalDiff = 0;

    for (int i = 0; i < k; i++)
    {
        actInstance.clear();
        calcInstances(actInstance, X, i, S);
        calcCentroidCoords(actCoords, actInstance, X);
        actDiff = calcIntraDiff(actCoords, actInstance, X);
        totalDiff += actDiff;
    }

    return totalDiff / k;
}

/**
 * @brief Get the fitness of the solution
 * @param S         Vector that contains the a solution
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param lambda    Value to be used for computing the fitness value.
 * @param k         Number of Centroids.
 * @return Value of the fitness obtained  
 */
double getFitness(vector<int> S, vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, double lambda, int k)
{
    return intraClusterDistance(S, X, k) + (lambda * infeasibility(S, ML, CL));
}

/**
 * @brief Update the current solution with the selected neighbor
 * @param S             Vector that contains the current solution.
 * @param change        "Virtual" Neighbor that makes S better than the current solution.
 * @param clusterList   Vector that mantains the number of instances assigned to each cluster.
 */
void updateSolution(vector<int> &S, pair<int,int> change, vector<int> &clusterList)
{
    clusterList[S[change.first]]--;
    clusterList[change.second]++;
    
    S[change.first] = change.second;
}

/**
 * @brief Use a Local Search Algorithm to obtain a solution for the Clustering with Restrictions Problem
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 * @param lambda    Value to be used for computing the fitness value.
 */
vector<int> localSearch(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
{
    int instances = X.size(),
        iterations = 0;

    // Vector that holds the  i instances associated with their cluster.
    vector<int> S (instances, -1);

    // Auxiliar vectors, to check if the solution has changed and to apply a neighbor to the current one.
    vector<int> oldS, actS;

    double bestFitness,
           actFitness;

    // Vector that holds the virtual neighborhood has a pair [l, i] that would make S[l]=i
    vector<pair<int,int>> virtuaNeighbors;
    pair<int, int> actualNeighbor;

    // This vector mantains the count of how many instances are assigned to a cluster.
    vector<int> clusterList (k);

    // Generate the initial solution and initial fitness.
    genInitSolution(S, k, clusterList);
    bestFitness = getFitness(S, X, ML, CL, lambda, k);

    do
    {
        oldS = S;

        // Generate the virtual neighborhood.
        virtuaNeighbors = genNeighbor(S, k, clusterList);

        // Iterate over the neighborhood.
        for (int i = 0; i < virtuaNeighbors.size(); i++)
        {

            // Pick one neighbor
            actualNeighbor = virtuaNeighbors[i];

            actS = S;
            actS[actualNeighbor.first] = actualNeighbor.second;

            // Calculate the fitness of changing the solution with that neighbor.
            actFitness = getFitness(actS, X, ML, CL, lambda, k);
            iterations++;

            // If its better, replace the solution with the change.
            if(actFitness < bestFitness)
            {
                updateSolution(S, actualNeighbor, clusterList);
                bestFitness = actFitness;
                break;
            }

            if(iterations >= 100000)
            {
                break;
            }
        }
                        
    }while (iterations < 100000 && S != oldS);
    
    return S;
}

/**
 * @brief Main function
 * @param n             Number of instances
 * @param d             Number of dimensions
 * @param k             Number of clusters
 * @param setPath       Path to the '*.dat' file
 * @param constPath     Path to the '*.const' file
 * @param randomSeed    Intenger to be used as the seed for the RNG
 * @return Execution status
 */
int main(int argc, char * argv[])
{

    if(argc != 7)
    {
        cout<<"\nThis program needs 6 parameters, it received "<<argc-1<<"."<<endl;
        cout<<"\nUsage: ./BL_exe n d k setPath constPath randomSeed"<<endl;
        cout<<"n            Number of instances."<<endl;
        cout<<"d            Dimensions of instances."<<endl;
        cout<<"k            Number of clusters."<<endl;
        cout<<"setPath      Path to the '*.dat' file"<<endl;
        cout<<"constPath    Path to the '*.const' file"<<endl;
        cout<<"randomSeed   Intenger to be used as the seed for the RNG"<<endl;
        exit(-1);
    }

    //  Declare the "stopwatch"
    clock_t timeBefore, timeAfter;

    // Getting data from input
    int numberInstances = stoi(argv[1]),
        dimensions = stoi(argv[2]),
        numClusters = stoi(argv[3]);

    string setPath = argv[4],
            constPath = argv[5];

    // Initializing the seed
    Set_random(stoi(argv[6]));


    // Defining main data structures
    vector<vector<float>> X;
        X.resize(numberInstances);
        for(int i = 0; i < numberInstances; i++) X[i].resize(dimensions);

    vector<vector<int>> MR;
        MR.resize(numberInstances);
        for(int i=0; i<numberInstances;i++) MR[i].resize(numberInstances);
    
    vector<vector<float>> centroids;
        centroids.resize(numClusters);
        for(int i = 0; i < numClusters; i++) centroids[i].resize(dimensions);

    vector<triplet> LR;
    vector<triplet> ML, CL;


    // Reading the files and filling up the data structures.
    loadData(setPath, constPath, X, MR);
    fillTriplet(numberInstances, MR, LR, ML, CL);

    // Generate the Lambda Constant.
    double lambda = genLambda(numberInstances, dimensions, X, LR.size());

    timeBefore = clock();
    
    //  Use a Local Search approach to the problem
    vector<int> result = localSearch(X, ML, CL, numClusters, lambda);

    timeAfter = clock();


    // ANALYSIS
    //  Gather the data
    double time = (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC;
    double intraClusterDist = intraClusterDistance(result, X, numClusters);
    int infeas = infeasibility(result, ML, CL);
    double fitness = intraClusterDist + (lambda * infeas);

    // Print the data
    cout<<"Fitness: "<< fitness <<"\tInfeas: "<<infeas<<"\tK Dist: "<<intraClusterDist<<"\tTime: "<<time<<"\tSeed: "<<stoi(argv[6])<<endl;


    return 0;
}

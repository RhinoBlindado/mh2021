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
     * Lambda = ceiling(D) / |R|
     *  - D is the largest distance between two points. The actual value has to be the ceiling of the distance.
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
 * @brief Fill the centroids dimensions with a random float.
 * @param centroids Array of k length with d dimensions.
 */
void fillCentroids(vector<vector<float>> &centroids)
{
    int k = centroids.size(), d = centroids[0].size();

    for(int i = 0; i < k; i++)
    {
        for(int j = 0; j < d; j++)
        {
            centroids[i][j] = Rand();
        }
    }
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
 * @brief Calculates the infeasibility of the current solution using triplets.
 * @param S     Vector of the current solution.
 * @param ML    Vector of MUST-LINK Triplets.
 * @param CL    Vector of CANNOT-LINK Triplets.
 * @return  The count of how many times a restrictions has been violated.
 */
int infeasibilityBL(vector<int> S, vector<triplet> ML, vector<triplet> CL)
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
 * @brief Random Wrapper for the RNG Library to be used in the Random Shuffle function
 * @param i     Index that is used for Random Shuffle to generate a number. 
 */
int rWrapper(int i)
{
    return Randint(0, i);
}

/**
 * @brief Check if a given instance in a centroid has a ML or CL restriction with other member of that centroid
 * @param MR    The Restriction Matrix
 * @param C     The vector that holds what instances are assigned to what cluster.
 * @param x_i   The ith component of the Restriction Matrix
 * @param x_j   The jth component of the Restriction Matrix
 * @param actC  The current cluster to be checked
 * @return 1 if a restriction is violated, 0 otherwise.
 */
int V(vector<vector<int>> MR, vector<int> C, int x_i, int x_j, int actC)
{
    bool checkR = false;

    // Check Must-Link
    if (MR[x_i][x_j] == 1 && C[x_j] != actC)
    {
        checkR = true;
    }

    // Check Cannot-Link
    if (MR[x_i][x_j] == -1 && C[x_j] == actC)
    {
        checkR = true;
    }
    return (int)checkR;
}

/**
 * @brief Check the infeasibility of asigning an instance to a particular cluster.
 * @param instance          Instance of X to be checked.
 * @param selectedCluster   A given cluster.
 * @param C                 The vector that holds what instances are assigned to what cluster.
 * @return The number of restrictions violated.
 */
int infeasibility(int instance, int selectedCluster, vector<int> C, vector<vector<int>> MR)
{
    int count = 0;
    int sizeC = C.size();

    for(int i; i < sizeC; i++)
    {
        if(C[i] == selectedCluster)
        {
            count += V(MR, C, instance, i, selectedCluster);
        }
    }

    return count;
}

/**
 * @brief Get the cluster with the minimal distance to a given instance
 * @param X                     Array that contains n instances with their d dimensions. 
 * @param centroids             The vector that holds all the centroid coordinates.           
 * @param selectedCentroids     The clusters with the lowest infeasibility.
 * @param instance              The current instance.
 * @return The cluster with the closest centroid to the instance
 */
int minDistance(vector<vector<float>> X, vector<vector<float>> centroids, vector<int> selectedCentroids, int instance)
{
    float maxDist = numeric_limits<float>::max();
    float actDist;
    int bestCentroid = -1;

    for (int i = 0; i < selectedCentroids.size(); i++)
    {
        actDist = 0;
        for(int j = 0; j < X[0].size(); j++)
        {
            actDist += pow(X[instance][j] - centroids[selectedCentroids[i]][j], 2);
        }

        actDist = sqrt(actDist);

        if(actDist < maxDist)
        {
            maxDist = actDist;
            bestCentroid = selectedCentroids[i];
        }
    }
    return bestCentroid;
}

/**
 * @brief Updates the coordinates of a given centroid.
 * @param actK      The current cluster
 * @param centroids The vector that holds all the centroid coordinates
 * @param X         Array that contains n instances with their d dimensions.
 * @param C_list    The vector that holds what instances are assigned to what cluster.       
 */
void updateDistance(int actK, vector<vector<float>> &centroids, vector<vector<float>> X, vector<int> C_list)
{
    vector<int> actList;
    int size = C_list.size();

    // Obtain the instances assigned to the cluster actK
    for (int i = 0; i < size; i++)
    {
        if(C_list[i] == actK)
        {
            actList.push_back(i);
        }
    }

    int dimensions = X[0].size();
    int instances = actList.size();

    float sum;
    
    // Obtain the average distance for each dimension
    for (int i = 0; i < dimensions; i++)
    {
        sum = 0;
        for(int j = 0; j < instances; j++)
        {
            sum += X[actList[j]][i];
        }

        if(instances > 0)
            centroids[actK][i] = sum / instances;
    }

}

/**
 * @brief Use a Greedy Algorithm to obtain a solution for the Clustering with Restrictions Problem
 * @param X         Array that contains n instances with their d dimensions.
 * @param MR        Array of restrictions
 * @param k         Number of clusters
 * @param centroids The d coordinates of the centroids of each k cluster.
 * @returns Vector of instances assinged to each cluster, ith instance is asigned to the C[i] cluster.
 */
vector<int> greed(vector<vector<float>> X, vector<vector<int>> MR, int k, vector<vector<float>> centroids)
{
    int instances = X.size(), maxInfeas;

    // This vector holds what instances are assigned to what cluster, begins with all -1s.
    vector<int> C (instances, -1);

    // This vector holds the infeasibility of each cluster.
    vector<int> C_infe (k, -1);

    // Copy of the C from the last iteration to check if anything has changed.
    vector<int> C_old;

    // This vector holds the clusters with the lowest infeasibility.
    vector<int> C_lowest;

    // Indexes of the X matrix.
    vector<int> RSI(instances);
        for (int i = 0; i < instances; i++) RSI[i] = i;

    //  Shuffle the indexes of the X matrix
    random_shuffle(RSI.begin(), RSI.end(), rWrapper);

    do
    {
        C_old = C;

        // Check the infeasibility of asigning i to the cluster j; that infeasibility is saved in C_infe[j]
        for(int i=0; i<instances; i++)
        {
            for(int j=0; j<k; j++)
            {
                C_infe[j] = infeasibility(RSI[i], j, C, MR);
            }

            // Assign the cluster with the less infeasibility or if equal score, the closest.
            maxInfeas = numeric_limits<int>::max();
            for (int j = 0; j < k; j++)
            {
                if(C_infe[j] <= maxInfeas)
                {
                    if(maxInfeas != C_infe[j])
                    {
                        C_lowest.clear();
                    }
                    C_lowest.push_back(j);
                    maxInfeas = C_infe[j];
                }
            }
            
            // If there is more than one cluster with the lowest infeasibility, check which one is the closest.
            if(C_lowest.size() > 1)
            {
                int aux = minDistance(X, centroids, C_lowest, RSI[i]);
                C[RSI[i]] = aux;
            }
            else
            {
                C[RSI[i]] = C_lowest.front();
            }
            C_lowest.clear();            
        }
            
        // Update centroids coordinates.
        for (int i = 0; i < k; i++)
        {
            updateDistance(i, centroids, X, C);
        }

    }while (C_old != C);

    return C;
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

    // Setup the centroids
    fillCentroids(centroids);  

    timeBefore = clock();
    
    //  Use a greedy approach to the problem
    vector<int> result = greed(X, MR, numClusters, centroids);

    timeAfter = clock();

    // ANALYSIS

    //  Gather the data
    double time = (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC;
    double intraClusterDist = intraClusterDistance(result, X, numClusters);
    double lambda = genLambda(numberInstances, dimensions, X, LR.size());
    int infeas = infeasibilityBL(result, ML, CL);
    double fitness = intraClusterDist + (lambda * infeas);

    // Print the data
    cout<<"Fitness: "<< fitness <<"\tInfeas: "<<infeas<<"\tK Dist: "<<intraClusterDist<<"\tTime: "<<time<<"\tSeed: "<<stoi(argv[6])<<endl;
    return 0;
}

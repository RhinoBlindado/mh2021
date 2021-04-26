// INTRODUCTION

/**
 *  [CASTELLANO]
 * 
 * Practica 2.b - Técnicas de Búsqueda basadas en Poblaciones para el Problema del Agrupamiento con Restricciones
 * Asignatura: Metaheuristicas
 * Autor: Valentino Lugli (Github: @RhinoBlindado)
 * Abril, Mayo 2021
 */

/**
 *  [ENGLISH]
 *
 * Practice 2.b - Population-based search techniques for the Clustering Problem with Restrictions
 * Course: Metaheuristics
 * Author: Valentino Lugli (Github: @RhinoBlindado)
 * April, May 2021
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

void generateInitialPop(vector<vector<int>> &pop, int k)
{
    for (int i = 0; i < pop.size(); i++)
    {
        for (int j = 0; j < k; j++)
        {
            pop[i][j] = j;
        }

        for(int j = k; j < pop[0].size(); j++)
        {
            pop[i][j] = Randint(0, k-1);
        }

    }

    for (int i = 0; i < pop.size(); i++)
    {
        random_shuffle(pop[i].begin(), pop[i].end(), rWrapper);
    }

}

void evaluateInitialPop()
{
    
}


/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
void AGG_UN(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
{

    vector<vector<int>> population;
        population.resize(50);
        for(int i = 0; i < 50; i++) population[i].resize(X.size());


    generateInitialPop(population, k);
    evaluateInitialPop();

/*    for (int i = 0; i < 100000; i++)
    {
       
    }
  */
}

/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
void AGG_SF(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
{
       
}

/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
void AGE_UN(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
{
       
}

/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
void AGE_SF(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
{
       
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
    //clock_t timeBefore, timeAfter;

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

    vector<triplet> LR;
    vector<triplet> ML, CL;


    // Reading the files and filling up the data structures.
    loadData(setPath, constPath, X, MR);
    fillTriplet(numberInstances, MR, LR, ML, CL);

    // Generate the Lambda Constant.
    double lambda = genLambda(numberInstances, dimensions, X, LR.size());

    
    /**
     *      GENETIC ALGORITHMS
     * 
     */ 
 //   timeBefore = clock();

    AGG_UN(X, ML, CL, numClusters, lambda);

 //   timeAfter = clock();


    // ANALYSIS
    //  Gather the data
 /*   double time = (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC;
    double intraClusterDist = intraClusterDistance(result, X, numClusters);
    int infeas = infeasibility(result, ML, CL);
    double fitness = intraClusterDist + (lambda * infeas);

    // Print the data
    cout<<"Fitness: "<< fitness <<"\tInfeas: "<<infeas<<"\tK Dist: "<<intraClusterDist<<"\tTime: "<<time<<"\tSeed: "<<stoi(argv[6])<<endl;
*/

    return 0;
}

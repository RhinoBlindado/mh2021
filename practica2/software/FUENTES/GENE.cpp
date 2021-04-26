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

//  GLOBAL VARS
float   uniformCrossChance = 0.7,
        staticSegCrossChange = 1;

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

class cromosome
{
    public:
        vector<int> genes;
        float fitness;
        bool hasChanges;
        vector<int> clusters;

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

void generateInitialPop(vector<cromosome> &pop, int k)
{
    for (int i = 0; i < pop.size(); i++)
    {
        for (int j = 0; j < k; j++)
        {
            pop[i].genes[j] = j;
            pop[i].clusters[j]++;
        }

        for(int j = k; j < pop[0].genes.size(); j++)
        {
            pop[i].genes[j] = Randint(0, k-1);
            pop[i].clusters[pop[i].genes[j]]++;
        }

    }

    for (int i = 0; i < pop.size(); i++)
    {
        random_shuffle(pop[i].genes.begin(), pop[i].genes.end(), rWrapper);
    }

}

void evaluateInitialPop(vector<cromosome> &pop, vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, double lambda, int k, pair<int, float> &bestCrom)
{
    for (int i = 0; i < pop.size(); i++)
    {
        pop[i].fitness = getFitness(pop[i].genes, X, ML, CL, lambda, k);
        
        if(pop[i].fitness < bestCrom.second)
        {
            bestCrom.first = i;
            bestCrom.second = pop[i].fitness;
        }
    }
    
}

void popSelection(vector<cromosome> pop, vector<cromosome> &parentPop)
{
    int fighterA,
        fighterB;

    for (int i = 0; i < pop.size(); i++)
    {
        fighterA = Randint(0, (pop.size()-1) );
        fighterB = Randint(0, (pop.size()-1) );

        if(pop[fighterA].fitness < pop[fighterB].fitness)
        {
            parentPop[i] = pop[fighterA];
        }
        else
        {
            parentPop[i] = pop[fighterB];
        }

    }
}

void popCrossUniform(vector<cromosome> &parentPop)
{
    vector<int> childA, 
                childB;

    int crossNum = (int)ceil(parentPop.size() * uniformCrossChance),
        geneChanges = (int)(parentPop[0].genes.size() / 2),
        geneSize = parentPop[0].genes.size() - 1,
        crossA,
        crossB;

    pair<int, int>  changeA, 
                    changeB;

    for (int i = 0; i < crossNum; i += 2)
    {
        childA = parentPop[i].genes;
        childB = childA;

        for (int j = 0; j < geneChanges; j++)
        {
            do
            {
                crossA = Randint(0, geneSize);
            }while(parentPop[i].clusters[ childA[crossA] ] > 1);

            do
            {
                crossB = Randint(0, geneSize);
            }while(parentPop[i].clusters[ childB[crossB] ] > 1);

            changeA.first = childA[crossA];
            changeA.second = parentPop[i+1].genes[crossA];

            changeB.first = childB[crossB];
            changeB.second = parentPop[i+1].genes[crossB];

            childA[crossA] = parentPop[i+1].genes[crossA];
            childB[crossB] = parentPop[i+1].genes[crossB];
        }

        parentPop[i].genes = childA;
        parentPop[i].hasChanges = true;
        parentPop[i].clusters[ changeA.first ]--;
        parentPop[i].clusters[ changeA.second ]++;

        parentPop[i+1].genes = childB;
        parentPop[i+1].hasChanges = true;
        parentPop[i+1].clusters[ changeB.first ]--;
        parentPop[i+1].clusters[ changeB.second ]++;
    }
}

void popMutation(vector<cromosome> &parentPop)
{
    int tVirus = Randint(0, parentPop[0].genes[0].size() - 1 ), 
        wesker = Randint(0, parentPop.size() - 1);

    
}


/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
void AGG_UN(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda, int popSize, bool operatorType)
{

    /*vector<vector<int>> population, parentPop;
        population.resize(popSize);
        for(int i = 0; i < popSize; i++) population[i].resize(X.size());  
        parentPop.resize(popSize);

    vector<float> popFitness (popSize, -1);
    vector<bool> needEval (popSize, true);*/

    vector<cromosome> population, parentPop, childPop;
    population.resize(popSize);
    parentPop.resize(popSize);
    
    for(int i = 0; i < popSize; i++) 
    {
        population[i].genes.resize(X.size());
        population[i].fitness = -1;
        population[i].hasChanges = false;
        population[i].clusters.resize(k);
    }

    pair<int, float> bestCrom(-1, numeric_limits<float>::max());

    generateInitialPop(population, k);
    evaluateInitialPop(population, X, ML, CL, lambda, k, bestCrom);

//    for (int i = 0; i < 100000; i++)
//    {
        popSelection(population, parentPop);
        popCrossUniform(parentPop);
        popMutation(parentPop);
//    }
}

/*        cout<<"OG POP:"<<endl;
        for (int i = 0; i < population.size(); i++)
        {
            for (int j = 0; j < population[0].genes.size(); j++)
            {
                cout<<population[i].genes[j]<<" ";
            }
            cout<<endl;
            
        }
*/

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
        numClusters = stoi(argv[3]),
        populationSize = 50;

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

    AGG_UN(X, ML, CL, numClusters, lambda, populationSize, true);

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

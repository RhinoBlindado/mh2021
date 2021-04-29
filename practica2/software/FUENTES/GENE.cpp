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
//  Common used functions for the practices
#include "libs/commonFunctions.h"
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
const float uniformCrossChance = 0.7,
            staticSegCrossChance = 1,
            mutationChance = 0.001;

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
    private:
    vector<int> genes;
    float fitness;
    bool hasChanges;
    vector<int> clusters;

    public:

    //  Constructors
    cromosome()
    {
        this->genes(0);
        this->fitness = -1;
        this->hasChanges = true;
        this->clusters(0);
    }

    cromosome(vector<int> arg_genes, float arg_fitness, bool arg_hasChanges, vector<int> arg_clusters)
    {
        this->genes = arg_genes;
        this->fitness = arg_fitness;
        this->hasChanges = arg_hasChanges;
        this->clusters = arg_clusters;
    }

    // Destructor
    ~cromosome()
    {
        this->genes.clear();
        this->clusters.clear();
    }

    // Getters
    vector<int> getGenes()
    {
        return this->genes;
    }

    int getGeneSize()
    {
        return this->genes.size();
    }

    int getGeneCluster(const int arg_gene)
    {
        return this->genes[arg_gene];
    }

    float getFitness()
    {
        return this->fitness;
    }

    bool getChanges()
    {
        return this->hasChanges;
    }

    vector<int> getClusters()
    {
        return this->clusters;
    }

    // Setters
    void initGenes(const int arg_geneSize)
    {
        this->genes.clear();
        this->genes.resize(arg_geneSize);
    }

    void initClusters(const int arg_clusterSize)
    {
        this->clusters.clear();
        this->clusters.resize(arg_clusterSize);
    }

    void setGenes(const vector<int> arg_genes)
    {
        this->genes.clear();
        this->genes = arg_genes;
    }

    void setFitness(const float arg_fitness)
    {
        this->fitness = arg_fitness;
    }

    void setChanges(const bool arg_hasChanges)
    {
        this->hasChanges = arg_hasChanges;
    }

    void setClusters(const vector<int> arg_clusters)
    {
        this->clusters = arg_clusters;
    }

    void changeClusterCount(const int arg_cluster, const int arg_count)
    {
        this->clusters[arg_cluster] += arg_count;
    }

    void changeGene(const int arg_index, const int arg_cluster)
    {
        this->genes[arg_index] = arg_cluster;
    }


    // Auxiliar Functions
    void printContents()
    {
        cout<<"---/Cromosome/---"<<endl;
        cout<<"\tGenes: [";
        for (int i = 0; i < genes.size(); i++)
        {
            cout<<genes[i]<<" ";
        }
        cout<<"]"<<endl;

        cout<<"\tFitness: "<<fitness<<endl;
        cout<<"\tHas changes? "<<hasChanges<<endl;
        cout<<"\tCluster count: [";

        for (int i = 0; i < clusters.size(); i++)
        {
            cout<<clusters[i]<<" ";
        }
        cout<<"]"<<endl;
        cout<<"---/Cromosome End/---"<<endl;
    }

    void printGenes()
    {
        cout<<"Genes: [ ";
        for (int i = 0; i < genes.size(); i++)
        {
           cout<<"["<<i<<"]="<<genes[i]<<"\t";
        }
        cout<<"]"<<endl;
        
    }

    void printClusters()
    {
        cout<<"Clusters: [ ";
        for (int i = 0; i < clusters.size(); i++)
        {
           cout<<"["<<i<<"]="<<clusters[i]<<"\t";
        }
        cout<<"]"<<endl;
    }

};

struct bestCromosome
{
    int pos;
    float fitness;
    bool isAlive;
    cromosome data;
};

// FUNCTIONS


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
//    cout<<"Infeas begin..."<<endl;

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
//    cout<<"Infeas done..."<<endl;
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
//     cout<<"Intracluster begin..."<<endl;
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
//    cout<<"Intracluster done..."<<endl;
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
    float iCD, infeas;
    
    iCD = intraClusterDistance(S, X, k);
    infeas = (lambda * infeasibility(S, ML, CL));

    return  iCD+infeas;
}

void generateInitialPop(vector<cromosome> &pop, int k)
{
    vector<int> auxGenes, auxClusters;
    auxGenes.resize(pop[0].getGeneSize());
    auxClusters.resize(k);

    for (int i = 0; i < pop.size(); i++)
    {

        for (int j = 0; j < k; j++)
        {
            auxClusters[j] = 0;

            auxGenes[j] = j;
            auxClusters[j]++;
        }

        for(int j = k; j < pop[0].genes.size(); j++)
        {
            auxGenes[j] = Randint(0, k-1);
            auxClusters[auxGenes[j]]++;
        }

        random_shuffle(auxGenes.begin(), auxGenes.end(), rWrapper);

        pop[i].setGenes(auxGenes);
        pop[i].setClusters(auxClusters);
    }
}

int evaluatePop(vector<cromosome> &pop, vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, double lambda, int k, bestCromosome &bestCrom)
{
    int internalEvaluations = 0;
    bestCrom.pos = -1;
    bestCrom.fitness = numeric_limits<float>::max();
    bestCrom.isAlive = true;

    for (int i = 0; i < pop.size(); i++)
    {
        if (pop[i].getChanges())
        {
            pop[i].setFitness( getFitness(pop[i].getGenes(), X, ML, CL, lambda, k) );
            pop[i].setChanges(false) ;
            internalEvaluations++;

        }
        
        if(pop[i].getFitness() < bestCrom.fitness)
        {
            bestCrom.pos = i;
            bestCrom.fitness = pop[i].getFitness();
        }
    }

    bestCrom.data = pop[ bestCrom.pos ];

    return internalEvaluations;
}

void popSelection(const vector<cromosome> pop, vector<cromosome> &parentPop, bestCromosome &bestCrom)
{
    int fighterA,
        fighterB,
        popSize = pop.size() - 1;

    bestCrom.isAlive = false;

    for (int i = 0; i < pop.size(); i++)
    {
        fighterA = Randint(0, popSize);
        fighterB = Randint(0, popSize);


        if(pop[fighterA].getFitness() < pop[fighterB].getFitness())
        {
            parentPop[i].setGenes( pop[fighterA].getGenes() );

            if(fighterA == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }
        else
        {
            parentPop[i].setGenes( pop[fighterB].getGenes() );

            if(fighterB == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }

    }
}

vector<cromosome> popCrossUniform(const vector<cromosome> parentPop, bestCromosome &bestCrom)
{
    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * uniformCrossChance),
        geneChanges = (int)(parentPop[0].getGeneSize() / 2),
        geneSize = parentPop[0].getGeneSize() - 1,
        crossA,
        crossB;

    if (bestCrom.isAlive && bestCrom.pos < crossNum)
    {
        bestCrom.isAlive = false;
    }
    
    for (int i = 0; i < crossNum; i += 2)
    {
        childA = parentPop[i].getGenes();
        clustA = parentPop[i].getClusters();

        childB = parentPop[i+1].getGenes();
        clustB = parentPop[i+1].getClusters();

        for (int j = 0; j < geneChanges; j++)
        {     
            do
            {
                crossA = Randint(0, geneSize);
            }while(clustA[ childA[crossA] ] == 1);
            do
            {
                crossB = Randint(0, geneSize);
            }while(clustB[ childB[crossB] ] == 1);

            clustA[ childA[crossA] ]--;
            clustA[ parentPop[i+1].genes[crossA] ]++;

            clustB[ childB[crossB] ]--;
            clustB[ parentPop[i].genes[crossB] ]++;

            childA[crossA] = parentPop[i+1].getGeneCluster(crossA);
            childB[crossB] = parentPop[i].getGeneCluster(crossB);
        }

        childPop[i].setGenes(childA);
        childPop[i].setChanges(true);
        childPop[i].setClusters(clustA);
        childPop[i].setFitness(-1);

        childPop[i+1].setGenes(childB);
        childPop[i+1].setChanges(true);
        childPop[i+1].setClusters(clustB);
        childPop[i+1].setFitness(-1);
    }

    return childPop;
}

vector<cromosome> popCrossFixedSeg(const vector<cromosome> parentPop, bestCromosome &bestCrom)
{
/*    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * uniformCrossChance),
        geneChanges = (int)(parentPop[0].genes.size() / 2),
        geneSize = parentPop[0].genes.size() - 1,
        crossA,
        crossB;
    
    if (bestCrom.isAlive && bestCrom.pos < crossNum)
    {
        bestCrom.isAlive = false;
    }
    
    for (int i = 0; i < crossNum; i += 2)
    {
        
    }

*/
}

void popMutation(vector<cromosome> &parentPop, bestCromosome bestCrom)
{
    int mutatedCroms = (int)ceil(mutationChance * parentPop.size() * parentPop[0].getGeneSize()),
        wesker = Randint(0, (parentPop.size() - mutatedCroms - 1)),
        clusterSize = parentPop[0].clusters.size() - 1,
        geneSize = parentPop[0].genes.size() - 1,
        mutationPos,
        tVirus;


    if(bestCrom.isAlive && (bestCrom.pos >= wesker && bestCrom.pos < wesker + mutatedCroms))
    {
        bestCrom.isAlive = false;
    }
    
    for (int i = 0; i < mutatedCroms; i++)
    {
        do
        {
            mutationPos = Randint(0, geneSize);
        }while(parentPop[i + wesker].clusters[ parentPop[i + wesker].genes[mutationPos] ] == 1 );

        do
        {
            tVirus = Randint(0, clusterSize);
        }while(parentPop[i + wesker].genes[mutationPos] == tVirus);

        parentPop[i + wesker].hasChanges = true;
        parentPop[i + wesker].fitness = -1;
        parentPop[i + wesker].clusters[tVirus]++;
        parentPop[i + wesker].clusters[ parentPop[i + wesker].genes[mutationPos] ]--;
        parentPop[i + wesker].genes[mutationPos] = tVirus;
    }
    
}

void replacePop(vector<cromosome> &childPop, vector<cromosome> parentPop, bestCromosome bestCrom)
{

    for (int i = 0; i < parentPop.size(); i++)
    {
        childPop[i] = parentPop[i];        
    }

    if (!bestCrom.isAlive)
    {
        childPop[0] = bestCrom.data;
    }
        
}

/**
 * @brief 
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.    
 */
vector<int> AGG(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda, int popSize, bool useUniformOperator)
{
    int evaluations = 0;

    vector<cromosome> population, parentPop, childPop;
    population.resize(popSize);
    parentPop.resize(popSize);
    
    for(int i = 0; i < popSize; i++) 
    {
        population[i].initGenes(X[0].size());
        population[i].initClusters(k);
        population[i].setFitness(-1);
        population[i].setChanges(true);
    }

    bestCromosome bestCrom;
    generateInitialPop(population, k);
    evaluatePop(population, X, ML, CL, lambda, k, bestCrom);

    while(evaluations < 100000)
    {
        popSelection(population, parentPop, bestCrom);

        //cout<<"Pop Cross"<<endl;
        if(useUniformOperator)
        {
           childPop = popCrossUniform(parentPop, bestCrom);
        }
        else
        {
            /*To-Do*/
            childPop = popCrossFixedSeg(parentPop, bestCrom);
        }

        popMutation(childPop, bestCrom);

        replacePop(population, childPop, bestCrom);

        evaluations += evaluatePop(population, X, ML, CL, lambda, k, bestCrom);

        cout<<"best="<<bestCrom.pos<<" Fitness="<<bestCrom.fitness<<endl;
   }

    return bestCrom.data.genes;

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
void AGE(vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int k, double lambda)
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
    cout<<argv[0]<<endl;
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
    timeBefore = clock();

    vector<int> result = AGG(X, ML, CL, numClusters, lambda, populationSize, true);

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

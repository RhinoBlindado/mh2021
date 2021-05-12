// INTRODUCTION

/**
 *  [CASTELLANO]
 * 
 * Practica 2.b - Técnicas de Búsqueda basadas en Poblaciones para el Problema del Agrupamiento con Restricciones
 * Asignatura: Metaheuristicas
 * Autor: Valentino Lugli (Github: @RhinoBlindado)
 * Fecha: Abril, Mayo 2021
 */

/**
 *  [ENGLISH]
 *
 * Practice 2.b - Population-based search techniques for the Clustering Problem with Restrictions
 * Course: Metaheuristics
 * Author: Valentino Lugli (Github: @RhinoBlindado)
 * Fecha: April, May 2021
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
//  Randomness.
#include "libs/random.h"
//  List.
#include <list>

using namespace std;

//  GLOBAL VARS
const float generationalCrossChance = 0.7,
            stationaryCrossChance = 1,
            mutationChance = 0.1;

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
        this->fitness = -1;
        this->hasChanges = true;
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
    vector<int> getGenes() const
    {
        return this->genes;
    }

    int getGeneSize() const
    {
        return this->genes.size();
    }

    int getGeneCluster(const int arg_gene) const
    {
        return this->genes[arg_gene];
    }

    float getFitness() const
    {
        return this->fitness;
    }

    bool getChanges() const
    {
        return this->hasChanges;
    }

    vector<int> getClusters() const
    {
        return this->clusters;
    }

    int getClusterSize() const
    {
        return this->clusters.size();
    }

    int getClusterCount( int arg_index ) const
    {
        return this->clusters[arg_index];
    }

    // Setters
    void initGenes(int arg_geneSize)
    {
        this->genes.clear();
        this->genes.resize(arg_geneSize);
    }

    void initClusters(int arg_clusterSize)
    {
        this->clusters.clear();
        this->clusters.resize(arg_clusterSize);
    }

    void setGenes(vector<int> arg_genes)
    {
        this->genes.clear();
        this->genes = arg_genes;
    }

    void setFitness(float arg_fitness)
    {
        this->fitness = arg_fitness;
    }

    void setChanges(bool arg_hasChanges)
    {
        this->hasChanges = arg_hasChanges;
    }

    void setClusters(vector<int> arg_clusters)
    {
        this->clusters.clear();
        this->clusters = arg_clusters;
    }

    void changeClusterCount(int arg_cluster, int arg_count)
    {
        this->clusters[arg_cluster] += arg_count;
    }

    void changeGene(int arg_index, int arg_cluster)
    {
        this->fitness= -1;
        this->hasChanges = true;

        this->clusters[arg_cluster]++;
        this->clusters[ this->genes[arg_index] ]--;

        this->genes[arg_index] = arg_cluster;
    }


    // Auxiliar Functions
    void printContents() const
    {
        cout<<"---/Cromosome/---"<<endl;
        cout<<"\tGenes:  [";
        for (int i = 0; i < genes.size(); i++)
        {
            cout<<"i="<<i<<" ["<<genes[i]<<"] ";
        }
        cout<<"]"<<endl;

        cout<<"\tFitness: "<<fitness<<endl;
        cout<<"\tHas changes? "<<hasChanges<<endl;
        cout<<"\tCluster count: [";

        for (int i = 0; i < clusters.size(); i++)
        {
            cout<<"i="<<i<<" ["<<clusters[i]<<"] ";
        }
        cout<<"]"<<endl;
        cout<<"---/Cromosome End/---"<<endl;
    }

    void printGenes() const
    {
        cout<<"Genes: [ ";
        for (int i = 0; i < genes.size(); i++)
        {
           cout<<"["<<i<<"]="<<genes[i]<<"\t";
        }
        cout<<"]"<<endl;
        
    }

    void printClusters() const
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
    return (Randint(0, i))%i;
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
double getFitness(const vector<int> &S, const vector<vector<float>> &X, const vector<triplet> &ML, const vector<triplet> &CL, double lambda, int k)
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

        for(int j = k; j < pop[0].getGeneSize(); j++)
        {
            auxGenes[j] = (Randint(0, k))%k;
            auxClusters[auxGenes[j]]++;
        }

        random_shuffle(auxGenes.begin(), auxGenes.end(), rWrapper);
        
        pop[i].setGenes(auxGenes);
        pop[i].setClusters(auxClusters);
    }
}

int evaluatePopGen(vector<cromosome> &pop, vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, double lambda, int k, bestCromosome &bestCrom)
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
            pop[i].setChanges(false);
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

void popSelectionGen(const vector<cromosome> &pop, vector<cromosome> &parentPop, bestCromosome &bestCrom)
{
    int fighterA,
        fighterB,
        popSize = pop.size();

    bestCrom.isAlive = false;

    for (int i = 0; i < pop.size(); i++)
    {
        fighterA = (Randint(0, popSize)%popSize);
        fighterB = (Randint(0, popSize)%popSize);


        if(pop[fighterA].getFitness() < pop[fighterB].getFitness())
        {
            parentPop[i] = pop[fighterA];

            if(fighterA == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }
        else
        {
            parentPop[i] = pop[fighterB];

            if(fighterB == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }

    }
}

vector<cromosome> popCrossUniformGen(const vector<cromosome> &parentPop, bestCromosome &bestCrom)
{
    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * generationalCrossChance),
        geneChanges = (int)(parentPop[0].getGeneSize() / 2),
        geneSize = parentPop[0].getGeneSize(),
        crossA = 0,
        crossB = 0;

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
                crossA = (Randint(0, geneSize))%geneSize;
            }while(clustA[ childA[crossA] ] == 1);
            do
            {
                crossB = (Randint(0, geneSize))%geneSize;
            }while(clustB[ childB[crossB] ] == 1);

            clustA[ childA[crossA] ]--;
            clustA[ parentPop[i+1].getGeneCluster(crossA) ]++;

            clustB[ childB[crossB] ]--;
            clustB[ parentPop[i].getGeneCluster(crossB) ]++;

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

vector<cromosome> popCrossFixedSegGen(const vector<cromosome> &parentPop, bestCromosome &bestCrom)
{
    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * generationalCrossChance),
        geneSize = parentPop[0].getGeneSize(),
        segStart,
        segLen,
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

        segLen = ( Randint(0, geneSize) )%geneSize;
        segStart = ( Randint(0, geneSize) )%geneSize;

        for (int j = 0; j < segLen; j++)
        {
            crossA = Randint(0, geneSize)%2;
            crossB = Randint(0, geneSize)%2;

            if(crossA && clustA[ childA[ (j + segStart)%geneSize ] ] > 1)
            {
                clustA[ childA[ (j + segStart)%geneSize ] ]--;
                clustA[ parentPop[i+1].getGeneCluster( (j + segStart)%geneSize ) ]++;

                childA[ (j + segStart)%geneSize ] = parentPop[i+1].getGeneCluster( ( (j + segStart)%geneSize ) );
            }

            if(crossB && clustB[ childB[ (j + segStart)%geneSize ] ] > 1)
            {
                clustB[ childB[ (j + segStart)%geneSize ] ]--;
                clustB[ parentPop[i].getGeneCluster( (j + segStart)%geneSize ) ]++;

                childB[ (j + segStart)%geneSize ] = parentPop[i].getGeneCluster( ( (j + segStart)%geneSize ) );
            }

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

void popMutationGen(vector<cromosome> &parentPop, bestCromosome bestCrom)
{
    int mutatedCroms = (int)ceil((mutationChance / parentPop[0].getGeneSize()) * parentPop.size() * parentPop[0].getGeneSize()),
        wesker = (Randint(0, (parentPop.size() - mutatedCroms - 1)))%parentPop.size(),
        clusterSize = parentPop[0].getClusterSize() - 1,
        geneSize = parentPop[0].getGeneSize(),
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
            mutationPos = (Randint(0, geneSize))%geneSize;
        }while(parentPop[i + wesker].getClusterCount( parentPop[i + wesker].getGeneCluster(mutationPos) ) == 1 );

        do
        {
            tVirus = (Randint(0, clusterSize))%clusterSize;
        }while(parentPop[i + wesker].getGeneCluster(mutationPos) == tVirus);

        parentPop[i + wesker].changeGene(mutationPos, tVirus);
    }
    
}

void replacePopGen(vector<cromosome> &childPop, vector<cromosome> parentPop, bestCromosome bestCrom)
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
vector<int> AGG(const vector<vector<float>> &X, const vector<triplet> &ML, const vector<triplet> &CL, int k, double lambda, int popSize, bool useUniformOperator)
{
    int evaluations = 0;

    vector<cromosome> population, parentPop, childPop;
    population.resize(popSize);
    parentPop.resize(popSize);
    
    for(int i = 0; i < popSize; i++) 
    {
        population[i].initGenes(X.size());
        population[i].initClusters(k);
        population[i].setFitness(-1);
        population[i].setChanges(true);
    }

    bestCromosome bestCrom;
    generateInitialPop(population, k);

    evaluatePopGen(population, X, ML, CL, lambda, k, bestCrom);
    while(evaluations < 100000)
    {
        popSelectionGen(population, parentPop, bestCrom);

        if(useUniformOperator)
        {
           childPop = popCrossUniformGen(parentPop, bestCrom);
        }
        else
        {
            childPop = popCrossFixedSegGen(parentPop, bestCrom);
        }

        popMutationGen(childPop, bestCrom);

        replacePopGen(population, childPop, bestCrom);
        evaluations += evaluatePopGen(population, X, ML, CL, lambda, k, bestCrom);
   }
    return bestCrom.data.getGenes();

}


int evaluatePopSta(vector<cromosome> &pop, const vector<vector<float>> &X, const vector<triplet> &ML, const vector<triplet> &CL, double lambda, int k, bestCromosome &bestCrom, vector<pair<int, float>> &worstTwoCroms)
{
    int internalEvaluations = 0,
        firstCrom = -1,
        secondCrom = -1;

    bool foundFirstWorst = false, 
         foundSecondWorst = false;

    float firstEval = -1,
          secondEval = -1;

    bestCrom.pos = -1;
    bestCrom.fitness = numeric_limits<float>::max();

    for (int i = 0; i < pop.size(); i++)
    {
        if (pop[i].getChanges())
        {
            pop[i].setFitness( getFitness(pop[i].getGenes(), X, ML, CL, lambda, k) );
            pop[i].setChanges(false);
            internalEvaluations++;
        }

        if (firstEval < pop[i].getFitness())
        {
            if (secondEval <= firstEval)
            {
                firstEval = secondEval;
                firstCrom = secondCrom;

                secondCrom = i;
                secondEval = pop[i].getFitness();
            }
            else
            {
                firstCrom = i;
                firstEval = pop[i].getFitness();
            }   
        }

        if(pop[i].getFitness() < bestCrom.fitness)
        {
            bestCrom.pos = i;
            bestCrom.fitness = pop[i].getFitness();
        }

    }

    bestCrom.data = pop[ bestCrom.pos ];

    worstTwoCroms[0].first = secondCrom;
    worstTwoCroms[0].second = secondEval;

    worstTwoCroms[1].first = firstCrom;
    worstTwoCroms[1].second = firstEval;

    return internalEvaluations;
}

void popSelectionSta(const vector<cromosome> pop, vector<cromosome> &parentPop)
{
    int fighterA,
        fighterB,
        popSize = pop.size();

    for (int i = 0; i < parentPop.size(); i++)
    {
        fighterA = (Randint(0, popSize)%popSize);
        fighterB = (Randint(0, popSize)%popSize);


        if(pop[fighterA].getFitness() < pop[fighterB].getFitness())
        {
            parentPop[i] = pop[fighterA];
        }
        else
        {
            parentPop[i] = pop[fighterB];
        }
    }
}

vector<cromosome> popCrossUniformSta(const vector<cromosome> &parentPop)
{
    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * stationaryCrossChance),
        geneChanges = (int)(parentPop[0].getGeneSize() / 2),
        geneSize = parentPop[0].getGeneSize(),
        crossA = 0,
        crossB = 0;

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
                crossA = (Randint(0, geneSize))%geneSize;
            }while(clustA[ childA[crossA] ] == 1);
            do
            {
                crossB = (Randint(0, geneSize))%geneSize;
            }while(clustB[ childB[crossB] ] == 1);

            clustA[ childA[crossA] ]--;
            clustA[ parentPop[i+1].getGeneCluster(crossA) ]++;

            clustB[ childB[crossB] ]--;
            clustB[ parentPop[i].getGeneCluster(crossB) ]++;

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

vector<cromosome> popCrossFixedSegSta(const vector<cromosome> &parentPop)
{
    vector<int> childA, clustA, 
                childB, clustB;

    vector<cromosome> childPop = parentPop;

    int crossNum = (int)ceil(parentPop.size() * stationaryCrossChance),
        geneSize = parentPop[0].getGeneSize(),
        segStart,
        segLen,
        crossA,
        crossB;
     
    for (int i = 0; i < crossNum; i += 2)
    {
        childA = parentPop[i].getGenes();
        clustA = parentPop[i].getClusters();

        childB = parentPop[i+1].getGenes();
        clustB = parentPop[i+1].getClusters();

        segLen = ( Randint(0, geneSize) )%geneSize;
        segStart = ( Randint(0, geneSize) )%geneSize;

        for (int j = 0; j < segLen; j++)
        {
            crossA = Randint(0, geneSize)%2;
            crossB = Randint(0, geneSize)%2;

            if(crossA && clustA[ childA[ (j + segStart)%geneSize ] ] > 1)
            {
                clustA[ childA[ (j + segStart)%geneSize ] ]--;
                clustA[ parentPop[i+1].getGeneCluster( (j + segStart)%geneSize ) ]++;

                childA[ (j + segStart)%geneSize ] = parentPop[i+1].getGeneCluster( ( (j + segStart)%geneSize ) );
            }

            if(crossB && clustB[ childB[ (j + segStart)%geneSize ] ] > 1)
            {
                clustB[ childB[ (j + segStart)%geneSize ] ]--;
                clustB[ parentPop[i].getGeneCluster( (j + segStart)%geneSize ) ]++;

                childB[ (j + segStart)%geneSize ] = parentPop[i].getGeneCluster( ( (j + segStart)%geneSize ) );
            }

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

void popMutationSta(vector<cromosome> &parentPop)
{
    int clusterSize = parentPop[0].getClusterSize() - 1,
        geneSize = parentPop[0].getGeneSize(),
        mutationPos,
        tVirus;
    
    for (int i = 0; i < parentPop.size(); i++)
    {
        
        if (fmod(Rand()+0.1, 1) > (1-mutationChance))
        {
            do
            {
                mutationPos = (Randint(0, geneSize))%geneSize;
            }while(parentPop[i].getClusterCount( parentPop[i].getGeneCluster(mutationPos) ) == 1 );

            do
            {
                tVirus = (Randint(0, clusterSize))%clusterSize;
            }while(parentPop[i].getGeneCluster(mutationPos) == tVirus);

            parentPop[i].changeGene(mutationPos, tVirus);
        }
    }
    
}

/**
 * @brief Auxiliar function that sorts the list with the pair that has the lowest fitness
 * @param a     First object to be compared.
 * @param b     Second object to be compared.
 * @returns TRUE if a has lower fitness than b, FALSE otherwise.
 */
bool compareFitness(const pair<int, float> &a, const pair<int, float> &b)
{
    return a.second < b.second;
}

/**
 * @brief Replacement Operator for the AGE
 * @param childPop  Vector of Cromosomes of the Child Population
 * @param parentPop Vector of Cromosomes of the Parent Population
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.  
 * @param lambda    Lambda value.           
 * @param worstTwoCroms Vector of pairs that cointains the two worst cromosomes of the generation.
 */
void replacePopSta(vector<cromosome> &childPop, vector<cromosome> parentPop, const vector<vector<float>> &X, const vector<triplet> &ML, const vector<triplet> &CL, double lambda, int k, vector<pair<int, float>> worstTwoCroms)
{
    list<pair<int, float>> genPop;
    pair<int, float> child1, child2;

    genPop.push_back(worstTwoCroms[0]);
    genPop.push_back(worstTwoCroms[1]);

    parentPop[0].setFitness( getFitness(parentPop[0].getGenes(), X, ML, CL, lambda, k) );
    parentPop[0].setChanges(false);

    parentPop[1].setFitness( getFitness(parentPop[1].getGenes(), X, ML, CL, lambda, k) );
    parentPop[1].setChanges(false);

    genPop.push_back( pair<int, float>(-1, parentPop[0].getFitness()) );
    genPop.push_back( pair<int, float>(-2, parentPop[1].getFitness()) );

    genPop.sort(compareFitness);

    child1 = genPop.front();
    genPop.pop_front();

    child2 = genPop.front();

    if (child1.first < 0)
    {
        childPop[ worstTwoCroms[1].first ] = parentPop[ abs(child1.first)-1 ]; 
    }

    if (child2.first < 0)
    {
        childPop[ worstTwoCroms[0].first ] = parentPop[ abs(child2.first)-1 ];
    }
}


/**
 * @brief           Genetic Algorithm with Stationary Population
 * @param X         Array that contains n instances with their d dimensions.
 * @param ML        Vector of MUST-LINK triplets.
 * @param CL        Vector of CANNOT-LINK triplets.
 * @param k         Number of Centroids.  
 * @param lambda    Lambda value
 * @param popSize   Size of Population
 * @param useUniformOperator    If TRUE, the Uniform Operator will be used, otherwise the Fixed Segment one will be used.  
 */
vector<int> AGE(const vector<vector<float>> &X, const vector<triplet> &ML, const vector<triplet> &CL, int k, double lambda, int popSize, bool useUniformOperator)
{
    int evaluations = 0;

    vector<cromosome> population, parentPop, childPop;
        population.resize(popSize);
        for(int i = 0; i < popSize; i++) 
        {
            population[i].initGenes(X.size());
            population[i].initClusters(k);
            population[i].setFitness(-1);
            population[i].setChanges(true);
        }
        parentPop.resize(2);

    bestCromosome bestCrom;

    vector<pair<int, float>> worstTwoCroms;
        worstTwoCroms.resize(2);

    generateInitialPop(population, k);
    evaluatePopSta(population, X, ML, CL, lambda, k, bestCrom, worstTwoCroms);

    while(evaluations < 100000)
    {
        popSelectionSta(population, parentPop);

        if(useUniformOperator)
        {
           childPop = popCrossUniformSta(parentPop);
        }
        else
        {
            childPop = popCrossFixedSegSta(parentPop);
        }

        popMutationSta(childPop);

        replacePopSta(population, childPop, X, ML, CL, lambda, k, worstTwoCroms);

        evaluatePopSta(population, X, ML, CL, lambda, k, bestCrom, worstTwoCroms);

        evaluations += 2;
   }
    
    return bestCrom.data.getGenes();
}

/**
 * @brief Print relevant information about a function.
 * @param   hint        A text string to be printed along with the data.
 * @param   S           Vector containing the resulting clustering.
 * @param   X           Matrix contaning n data with d dimensions.
 * @param   ML          Vector of MUST-LINK triplets.
 * @param   CL          Vector of CANNOT-LINK triplets.
 * @param   numClusters Number of clusters.
 * @param   lambda      Lambda value calculated.
 * @param   time        Time that takes to execute the function.
 * @param   seed        The seed used in the execution.
 */
void output(string hint, vector<int> S, vector<vector<float>> X, vector<triplet> ML, vector<triplet> CL, int numClusters, double lambda, double time, int seed)
{

    double intraClusterDist = intraClusterDistance(S, X, numClusters);
    int infeas = infeasibility(S, ML, CL);
    double fitness = intraClusterDist + (lambda * infeas);
    
    cout<<hint<<"\t"<<infeas <<"\t"<<intraClusterDist<<"\t"<<fitness<<"\t"<<time<<"\t"<<seed<<endl;
    
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
        numClusters = stoi(argv[3]),
        populationSize = 50;

    string setPath = argv[4],
            constPath = argv[5];

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

    // Evaluate the algorithms
    
    // Initializing the seed
    Set_random(stoi(argv[6]));
    timeBefore = clock();
    vector<int> AGG_UN = AGG(X, ML, CL, numClusters, lambda, populationSize, true);
    timeAfter = clock();

    output("AGG_UN", AGG_UN, X, ML, CL, numClusters, lambda, (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC, stoi(argv[6]));

    Set_random(stoi(argv[6]));
    timeBefore = clock();
    vector<int> AGG_SF = AGG(X, ML, CL, numClusters, lambda, populationSize, false);
    timeAfter = clock();

    output("AGG_SF", AGG_SF, X, ML, CL, numClusters, lambda, (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC, stoi(argv[6]));

    Set_random(stoi(argv[6]));
    timeBefore = clock();
    vector<int> AGE_UN = AGE(X, ML, CL, numClusters, lambda, populationSize, true);
    timeAfter = clock();

    output("AGE_UN", AGE_UN, X, ML, CL, numClusters, lambda, (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC, stoi(argv[6]));

    Set_random(stoi(argv[6]));
    timeBefore = clock();
    vector<int> AGE_SF = AGE(X, ML, CL, numClusters, lambda, populationSize, false);
    timeAfter = clock();

    output("AGE_SF", AGE_SF, X, ML, CL, numClusters, lambda, (double)(timeAfter-timeBefore) / CLOCKS_PER_SEC, stoi(argv[6]));

    return 0;
}

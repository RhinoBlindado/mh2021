
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
//  Use UNIX time to get randomness.
#include <ctime>
//  Do the random shuffle.
#include <algorithm> 
//  Get max size of data types.
#include <limits>
//  Perform math calculations.
#include <cmath>
//  Take the time between executions.
#include <ctime>

using namespace std;

//  STRUCTS
struct triplet
{
    int x_0;
    int x_1;
    int R;
};

// FUNCTIONS
void loadData(string dataSetPath, string dataSetRestrictionsPath, vector<vector<float>> &X, vector<vector<int>> &MR)
{
    ifstream dataSetFile (dataSetPath);
    ifstream dataSetRestrictions (dataSetRestrictionsPath);
    string line, item;

//    float upBound, loBound;
    
    if(dataSetFile.is_open())
    {
        int i=0, j;
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
        int i=0, j;
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

float rWrapper(int i)
{
    return Randint(0, i);
}

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

int infeasibility(int instance, int selectedCentroid, vector<int> assignedList, vector<vector<int>> MR)
{
    int count = 0;
    int sizeCentroid = assignedList.size();

    for(int i; i<sizeCentroid; i++)
    {
        if(assignedList[i] == selectedCentroid)
        {
            count += V(MR, assignedList, instance, i, selectedCentroid);
        }
    }

    return count;
}

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

void updateDistance(int actK, vector<vector<float>> &centroids, vector<vector<float>> X, vector<int> C_list)
{
    vector<int> actList;
    int size = C_list.size();


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

vector<int> greed(vector<vector<float>> X, vector<vector<int>> MR, int k, vector<vector<float>> centroids)
{
    int instances = X.size();
    vector<int> RSI(instances);
        for (int i = 0; i < instances; i++) RSI[i] = i;

    vector<int> C (instances, -1);
    vector<int> C_infe (k, -1);
    vector<int> C_old;
    vector<int> C_lowest;

    int maxInfeas;

    //  Shuffle the indexes of the X matrix
    random_shuffle(RSI.begin(), RSI.end(), rWrapper);

    do
    {
        C_old = C;

        // Check the infeasibility of asigning i to the cluster j
        for(int i=0; i<instances; i++)
        {
            for(int j=0; j<k; j++)
            {
                C_infe[j] = infeasibility(RSI[i], j, C, MR);
            }

            // Assign the centroid with the less feasibility or if equal score, the closest.
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
            
            if(C_lowest.size() > 1)
            {
                int aux = minDistance(X, centroids, C_lowest, RSI[i]);
                C[RSI[i]] = aux;
            }
            else
            {
                C[RSI[i]] = C_lowest.front();
            }
        }
        C_lowest.clear();

        // Update centroids
        for (int i = 0; i < k; i++)
        {
            updateDistance(i, centroids, X, C);
        }

    } while (C_old != C);

    return C;
}


// 1 2 3 4       5                6
// n d k set.dat set_const.const  randomSeed
int main(int argc, char * argv[])
{
    // Initializing the seed
    Set_random(time(NULL));

    //  Declare the "stopwatch"
    clock_t timeBefore, timeAfter;

    // Getting data from input
    int numberInstances = stoi(argv[1]),
        dimensions = stoi(argv[2]),
        numClusters = stoi(argv[3]);

    string setPath = argv[4],
            constPath = argv[5];

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


    // Reading the files and filling up the data structures.
    loadData(setPath, constPath, X, MR);

    // Setup the centroids
    fillCentroids(centroids);  


    timeBefore = clock();
    
    //  Use a greedy approach to the problem
    vector<int> result = greed(X, MR, numClusters, centroids);

    timeAfter = clock();

    for (int i = 0; i < result.size(); i++)
    {
//        cout<<result[i]<<endl;
    }
    
    return 0;
}

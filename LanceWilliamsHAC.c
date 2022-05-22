// This program was written by ETHAN FONG
// on 16-04-2022

// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering


#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX

double LanceWilliamsFormula(double distIK, double distJK, int method);
void findMinimumDistance(double **distArray, int nV, int *I, int *J);
Dendrogram *initialiseDendrogramArray(Graph g);
Dendrogram newDendrogram(int v);
double **initialiseDistArray(Graph g);
void freeDistArray(double **distArray, int nV);
double vertexDistance(Graph g, Vertex v, Vertex w);

/*******************************************************************************
                                   FUNCTIONS
*******************************************************************************/

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
    int nV = GraphNumVertices(g);
    Dendrogram *dendArray = initialiseDendrogramArray(g);
    Dendrogram new;
    double **distArray = initialiseDistArray(g);
    for (int nClusters = nV; nClusters > 1; nClusters--) {
        // I and J will be modified to represent the two clusters to be merged
        int I = 0;
        int J = 0;
        findMinimumDistance(distArray, nV, &I, &J);
        // iterate through all remaining clusters k to update distance to new
        // cluster, using Lance William's Formula
        for (int k = 0; k < nV; k++) {
            // cells with value -1 will be ignored
            if (distArray[k][J] != -1 && distArray[k][J] != -1) {
                double distIK = distArray[I][k];
                double distJK = distArray[J][k];
                // integer J will be used to represent new cluster
                distArray[J][k] = LanceWilliamsFormula(distIK, distJK, method);
                distArray[k][J] = distArray[J][k];
            }
            // cluster I will be 'deleted'
            // all cells with cluster I will be flagged to be ignored
            distArray[k][I] = -1;
            distArray[I][k] = -1;
        }
        new = newDendrogram(-1);
        new->left = dendArray[J];
        new->right = dendArray[I];
        // new dendrogram will be represented by J
        // cell I will be 'deleted'
        dendArray[J] = new;
        dendArray[I] = NULL;
    }
    free(dendArray);
    freeDistArray(distArray, nV);
    // return last dendrogram created, i.e. dendrogram will all clusters
    return new;
}

/*******************************************************************************
                               HELPER FUNCTIONS
*******************************************************************************/

// Calculates the distance from new cluster IJ to cluster K using Lance-Williams 
// formula
// For single linkage and complete linkage cases, the formula simplifies to 
// min and max functions respectively
double LanceWilliamsFormula(double distIK, double distJK, int method) {
    double min = distIK;
    double max = distJK;
    if (distIK > distJK) {
        min = distJK;
        max = distIK;
    }
    if (method == SINGLE_LINKAGE || max == INFINITY) {
        // in the case the method is complete linkage and the max distance
        // between distIK and distJK is INFINITY, i.e. not connected by an edge
        // choose min
        // in the case that neither clusters I or J are connected to K, then
        // distance between IJ and K will be INFINITY
        return min;
    }
    return max;
}

// Calculates the minimum distance in the between two vertices in g by using the
// distance array
// Takes two integer pointers, I and J, and modifies them to the two vertices 
// that have the minimum distance
void findMinimumDistance(double **distArray, int nV, int *I, int *J) {
    double minDist = INFINITY;
    for (int i = 1; i < nV; i++) {
         // j is always less than i
        for (int j = 0; j < i; j++) {
            if (distArray[i][j] >= 0 && distArray[i][j] <= minDist) {
                // ignores dist values of -1
                // will eventually choose dist values of INFINITY,
                // i.e. unconnected clusters
                *I = i;
                *J = j;
                minDist = distArray[i][j];
            }
        }
    }
    return;
}

// Initialises a new Dendrogram array
// Returns a pointer to an array of dendrograms that are initialised to each 
// point to a single dendrogram leaf that corresponds to each vertex
Dendrogram *initialiseDendrogramArray(Graph g) {
    int nV = GraphNumVertices(g);
    Dendrogram *dendArray = calloc(nV, sizeof(*dendArray));
    if (dendArray == NULL) {
        fprintf(stderr, "error: couldn't allocate new dendrogram array\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nV; i++) {
        dendArray[i] = newDendrogram(i);
    }
    return dendArray;
}

// Creates a new Dendrogram
// Returns a pointer to a new dendrogram node of vertex v
Dendrogram newDendrogram(int v) {
    Dendrogram new = malloc(sizeof(*new));
    if (new == NULL) {
        fprintf(stderr, "error: couldn't allocate new dendrogram node\n");
        exit(EXIT_FAILURE);
    }
    new->vertex = v;
    new->left = NULL;
    new->right = NULL;
    return new;
}

// Initialises a new 2D distance array
// Returns a pointer to 2D array of doubles that represent the distance between
// the two vertices represented by the row and column indices.
double **initialiseDistArray(Graph g) {
    int nV = GraphNumVertices(g);
    double **new2DArray = calloc(nV, sizeof(*new2DArray)); 
    if (new2DArray == NULL) {
        fprintf(stderr, "error: couldn't allocate new distance array\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nV; i++) {
        new2DArray[i] = calloc(nV, sizeof(double));
        if (new2DArray[i] == NULL) {
            fprintf(stderr, "error: couldn't allocate new distance array\n");
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 0; i < nV; i++) {
        for (int j = 0; j < nV; j++) {
            if (i == j) {
                // the distances between the same vertex are initialised to -1
                new2DArray[i][j] = -1;
            } else {
                new2DArray[i][j] = vertexDistance(g, i, j);
            }
        }
    }
    return new2DArray;
}

// Frees all memory associated with the given 2D distance array
// The width and length of the array is nV
void freeDistArray(double **distArray, int nV) {
    for (int i = 0; i < nV; i++) {
        free(distArray[i]);
    }
    free(distArray);
}

// Calculates the distance between two vertices v and w
// Used to initialise the distance array, where each vertex is treated as 
// a cluster
double vertexDistance(Graph g, Vertex v, Vertex w) {
    double max_weight = 0;
    AdjList vOut = GraphOutIncident(g, v);
    // compare max edge weight to edges from v to w
    for (AdjList curr = vOut; curr != NULL; curr = curr->next) {
        if (curr->v == w && curr->weight > max_weight) {
            max_weight = curr->weight;
        }
    }
    AdjList wOut = GraphOutIncident(g, w);
    // compare max edge weight to edges from w to v
    for (AdjList curr = wOut; curr != NULL; curr = curr->next) {
        if (curr->v == v && curr->weight > max_weight) {
            max_weight = curr->weight;
        }
    }
    if (max_weight > 0) {
        // edge between v and w exists
        return 1 / max_weight;
    }
    // otherwise no edge exists
    return INFINITY;
}



// This program was written by ETHAN FONG
// on 16-04-2022

// Centrality Measures API implementation

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"

double sumShortestPaths(Graph g, Vertex v);
double numReachable(Graph g, Vertex v);
double numShortestPaths(Graph g, Vertex src, Vertex dest, NodeData *nodeArray);
double numVPaths(Graph g, Vertex v, Vertex src, Vertex dest, 
                  NodeData *nodeArray);
double betweenessCentralityOfV(Graph g, Vertex v);

/*******************************************************************************
                                   FUNCTIONS
*******************************************************************************/

double *closenessCentrality(Graph g) {
    int nV = GraphNumVertices(g);
    double *closenessArray = calloc(nV, sizeof(double));
    if (closenessArray == NULL) {
        fprintf(stderr, "error: couldn't allocate new closeness array\n");
        exit(EXIT_FAILURE);
    }
    // calculate closeness centrality for every vertex in g
    for (int i = 0; i < nV; i++) {
        double n = numReachable(g, i);
        double sumPaths = sumShortestPaths(g, i);
        if (nV == 1 || sumPaths == 0) {
            // in the case the denominator will become zero
            closenessArray[i] = 0;
        } else {
            closenessArray[i] = (n - 1) * (n - 1) / ((nV - 1) * sumPaths);
        }
    }
    return closenessArray;
}

double *betweennessCentrality(Graph g) {
    int nV = GraphNumVertices(g);
    double *betweenessArray = calloc(nV, sizeof(double));
    if (betweenessArray == NULL) {
        fprintf(stderr, "error: couldn't allocate new betweeness array\n");
        exit(EXIT_FAILURE);
    }
    // calculate betweeness centrality for every vertex in g
    for (int i = 0; i < nV; i++) {
        betweenessArray[i] = betweenessCentralityOfV(g, i);
    }
    return betweenessArray;
}

/*******************************************************************************
                               HELPER FUNCTIONS
*******************************************************************************/

// Calculates the sum of the shortest-path distances from vertex v to all other
// vertices
double sumShortestPaths(Graph g, Vertex v) {
    // create dijkstra array with all shortest-path distances
    NodeData *nodeArray = dijkstra(g, v);
    int nV = GraphNumVertices(g);
    double sum = 0;
    for (int i = 0; i < nV; i++) {
        if (i != v && nodeArray[i].dist < INFINITY) {
            sum += nodeArray[i].dist;
        }
    }
    freeNodeData(nodeArray, nV);
    return sum;
}

// Calculates the number of reachable vertices in Graph g from vertex v,
// including vertex v itself
double numReachable(Graph g, Vertex v) {
    NodeData *nodeArray = dijkstra(g, v);
    int nV = GraphNumVertices(g);
    double reachable = 0;
    for (int i = 0; i < nV; i++) {
        if (nodeArray[i].dist < INFINITY) {
            // vertex i is reachable
            reachable++;
        }
    }
    freeNodeData(nodeArray, nV);
    return reachable;
}

// Calculates the number of paths of shortest-path length from vertex src to 
// vertex dest in Graph g
// Takes a dijkstra array as argument nodeArray
//
// Recursively finds the number of paths to dest from src by calculating it as
// the sum of paths from src to dest vertex's predecessors
double numShortestPaths(Graph g, Vertex src, Vertex dest, NodeData *nodeArray) {
    if (src == dest) {
        return 1;
    }
    PredNode *curr = nodeArray[dest].pred;
    double sumPaths = 0;
    for (; curr != NULL; curr = curr->next) {
        // iterate through the predecessors of dest to sum the shortest paths
        // from src to them
        sumPaths += numShortestPaths(g, src, curr->v, nodeArray);
    }
    return sumPaths;
}

// Calculates the number of shortest paths from src to dest that contain 
// vertex v
double numVPaths(Graph g, Vertex v, Vertex src, Vertex dest, 
                  NodeData *nodeArray) {
    double fromV = numShortestPaths(g, v, dest, nodeArray);
    double toV = numShortestPaths(g, src, v, nodeArray);
    return fromV * toV;
}

// Calculates the betweeness centrality of a vertex v
double betweenessCentralityOfV(Graph g, Vertex v) {
    double bcOfV = 0;
    int nV = GraphNumVertices(g);
    // iterate through every pair of vertices i and j where i and j do not 
    // equal v
    for (int i = 0; i < nV; i++) {
        if (i == v) {
            continue;
        }
        NodeData *nodeArray = dijkstra(g, i);
        for (int j = 0; j < nV; j++) {
            if (j == v || j == i || nodeArray[j].dist == INFINITY) {
                // j equals i or v, or j is unreachable from i
                continue;
            }
            double nVPaths = numVPaths(g, v, i, j, nodeArray);
            double numPaths = numShortestPaths(g, i, j, nodeArray);
            bcOfV += nVPaths / numPaths;
        }
        freeNodeData(nodeArray, nV);
    }
    return bcOfV;
}





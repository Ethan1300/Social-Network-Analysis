// This program was written by ETHAN FONG
// on 14-04-2022

// Dijkstra API implementation

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

PredNode *PredNodeNew(Vertex v);
void PredNodeInsert(NodeData *data, Vertex pred);
void PredNodeFreeList (NodeData *data);

/*******************************************************************************
                                   FUNCTIONS
*******************************************************************************/


NodeData *dijkstra(Graph g, Vertex src) {
    int nV = GraphNumVertices(g);
    NodeData *nodeArray = calloc(nV, sizeof(NodeData));
    if (nodeArray == NULL) {
        fprintf(stderr, "error: couldn't allocate new dijkstra array\n");
        exit(EXIT_FAILURE);
    }
    PQ pq = PQNew();
    // Initialise NodeData array and priority queue
    for (int i = 0; i < nV; i++) {
        nodeArray[i].dist = INFINITY;	 
        nodeArray[i].pred = NULL;
        PQInsert(pq, i, INFINITY);
    }
    // Initialise source entry for NodeData array and priority queue
    nodeArray[src].dist = 0;
    PQUpdate(pq, src, 0);
    while (!PQIsEmpty(pq)) {
        // find vertex with minimum distance from source
        Vertex latest = PQDequeue(pq);
        if (nodeArray[latest].dist == INFINITY) {
            // no new vertices are reachable
            PQFree(pq);
            return nodeArray;
        }
        int predDist = nodeArray[latest].dist;
        AdjList newEdges = GraphOutIncident(g, latest);
        // loop through all outgoing edges from latest vertex and relax edges
        for (AdjList curr = newEdges; curr != NULL; curr = curr->next) {
            if (nodeArray[curr->v].dist > predDist + curr->weight) {
                // new path is shorter than shortest known path
                nodeArray[curr->v].dist = predDist + curr->weight;
                PredNodeFreeList(&nodeArray[curr->v]);
                PredNodeInsert(&nodeArray[curr->v], latest);
                PQUpdate(pq, curr->v, nodeArray[curr->v].dist);
            } else if (nodeArray[curr->v].dist == predDist + curr->weight) {
                // new path is same length as shortest known path
                PredNodeInsert(&nodeArray[curr->v], latest);
            }
        }
    }
    PQFree(pq);
    return nodeArray;
}

void freeNodeData(NodeData *data, int nV) {
    for (int i = 0; i < nV; i++) {
        PredNodeFreeList(&data[i]);
    }
    free(data);
}

/*******************************************************************************
                               HELPER FUNCTIONS
*******************************************************************************/

// Creates a new Prednode
//
//	Takes the vertex of the new PredNode
//	Returns a pointer to the new Prednode
PredNode *PredNodeNew(Vertex v) {
    PredNode *new = malloc(sizeof(*new));
    if (new == NULL) {
        fprintf(stderr, "error: couldn't allocate new PredNode\n");
        exit(EXIT_FAILURE);
    }
    new->v = v;
    new->next = NULL;
    return new;
}

// Inserts a PredNode into a list in ascending order
//
//	Takes a pointer to a NodeData struct that has the list and the vertex
//	being inserted
void PredNodeInsert(NodeData *data, Vertex pred) {
    PredNode *curr = data->pred;
    // list is empty
    if (curr == NULL) {
        data->pred = PredNodeNew(pred);
        return;
    }
    // new PredNode is the new head; new PredNode has smallest vertex
    if (pred < curr->v) {
        PredNode *new = PredNodeNew(pred);
        new->next = curr;
        data->pred = new;
        return;
    }
    // insert new PredNode after existing PredNode
    for (;curr->next != NULL && pred > curr->next->v; curr = curr->next) {
    }
    PredNode *new = PredNodeNew(pred);
    new->next = curr->next;
    curr->next = new;
    return;
}

// Frees all memory associated with the PredNode list of the given NodeData 
// struct
void PredNodeFreeList (NodeData *data) {
    while (data->pred != NULL) {
        PredNode *temp = data->pred;
        data->pred = temp->next;
        free(temp);
    }
    return;
}



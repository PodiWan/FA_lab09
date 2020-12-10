#include <cstdlib>
#include <cstring>
#include <queue>
#include <iostream>
#include "bfs.h"

int get_neighbors(const Grid *grid, Point p, Point neighb[]) {
    int iNeighbors = 0;

    if (p.col + 1 < grid->cols && grid->mat[p.row][p.col + 1] == 0)
        neighb[iNeighbors++] = Point(p.row, p.col + 1);

    if (p.row - 1 > 0 && grid->mat[p.row - 1][p.col] == 0)
        neighb[iNeighbors++] = Point(p.row - 1, p.col);

    if (p.col - 1 > 0 && grid->mat[p.row][p.col - 1] == 0)
        neighb[iNeighbors++] = Point(p.row, p.col - 1);

    if (p.row + 1 < grid->rows && grid->mat[p.row + 1][p.col] == 0)
        neighb[iNeighbors++] = Point(p.row + 1, p.col);

    return iNeighbors;
}

void grid_to_graph(const Grid *grid, Graph *graph) {
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node *nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (grid->mat[i][j] == 0) {
                nodes[i][j] = (Node *) malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            } else {
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node **) malloc(graph->nrNodes * sizeof(Node *));
    k = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (nodes[i][j] != NULL) {
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for (i = 0; i < graph->nrNodes; ++i) {
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if (graph->v[i]->adjSize != 0) {
            graph->v[i]->adj = (Node **) malloc(graph->v[i]->adjSize * sizeof(Node *));
            k = 0;
            for (j = 0; j < graph->v[i]->adjSize; ++j) {
                if (neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0) {
                    graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if (k < graph->v[i]->adjSize) {
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node **) realloc(graph->v[i]->adj, k * sizeof(Node *));
            }
        }
    }
}

void free_graph(Graph *graph) {
    if (graph->v != NULL) {
        for (int i = 0; i < graph->nrNodes; ++i) {
            if (graph->v[i] != NULL) {
                if (graph->v[i]->adj != NULL) {
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

void bfs(Graph *graph, Node *s, Operation *op) {
    // TODO: implement the BFS algorithm on the graph, starting from the node s
    // at the end of the algorithm, every node reachable from s should have the color BLACK
    // for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
    // for counting the number of operations, the optional op parameter is received
    // since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
    // if(op != NULL) op->count();

    if (op != NULL)
        op->count();

    for (int i = 0; i < graph->nrNodes; ++i) {
        graph->v[i]->color = COLOR_WHITE;
        graph->v[i]->dist = graph->nrNodes + 1;
        graph->v[i]->parent = nullptr;
        if (op != NULL)
            op->count(3);
    }
    s->color = COLOR_GRAY;
    s->dist = 0;
    s->parent = nullptr;
    if (op != NULL)
        op->count(3);

    std::queue<Node*> qNodes;
    qNodes.push(s);
    if (op != NULL)
        op->count();
    while (!qNodes.empty()) {
        Node* ptrAux = qNodes.front();
        if (op != NULL)
            op->count();
        qNodes.pop();
        for (int i = 0; i < ptrAux->adjSize; ++i) {
            if (op != NULL)
                op->count();
            if (ptrAux->adj[i]->color == COLOR_WHITE) {
                ptrAux->adj[i]->color = COLOR_GRAY;
                ptrAux->adj[i]->dist = ptrAux->dist + 1;
                ptrAux->adj[i]->parent = ptrAux;
                qNodes.push(ptrAux->adj[i]);
                if (op != NULL)
                    op->count(3);
            }
        }
        ptrAux->color = COLOR_BLACK;
        if (op != NULL)
            op->count();
    }
}

std::string prettyPrint(NodeR2* pR2Array, int iLevel){
    std::string strOut = "\n";

    for(auto i = 0; i < iLevel; ++i)
        strOut += "\t";
    strOut += "(" + std::to_string(pR2Array->position.row) + ", " + std::to_string(pR2Array->position.col) + ")";

    for (int i = 0; i < pR2Array->childrenSize; ++i)
        strOut += prettyPrint(pR2Array->children[i], iLevel + 1);

    return strOut;
}

void print_bfs_tree(Graph *graph) {
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int *p = NULL; //the parent array
    Point *repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int *transf = (int *) malloc(graph->nrNodes * sizeof(int));
    for (int i = 0; i < graph->nrNodes; ++i) {
        if (graph->v[i]->color == COLOR_BLACK) {
            transf[i] = n;
            ++n;
        } else {
            transf[i] = -1;
        }
    }
    if (n == 0) {
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int *) malloc(n * sizeof(int));
    repr = (Point *) malloc(n * sizeof(Node));
    for (int i = 0; i < graph->nrNodes && !err; ++i) {
        if (graph->v[i]->color == COLOR_BLACK) {
            if (transf[i] < 0 || transf[i] >= n) {
                err = 1;
            } else {
                repr[transf[i]] = graph->v[i]->position;
                if (graph->v[i]->parent == NULL) {
                    p[transf[i]] = -1;
                } else {
                    err = 1;
                    for (int j = 0; j < graph->nrNodes; ++j) {
                        if (graph->v[i]->parent == graph->v[j]) {
                            if (transf[j] >= 0 && transf[j] < n) {
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if (!err) {
        // TODO: pretty print the BFS tree
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
        auto pArray = new NodeR2*[n]();
        int* ptrCounts = new int[n]();

        for (int i = 0; i < n; ++i) {
            pArray[i] = (NodeR2*)malloc(sizeof(NodeR2));
            pArray[i]->position = repr[i];
            pArray[i]->childrenSize = 0;
        }

        for(int i = 0; i < n; ++i) {
            if(p[i] != -1)
                pArray[p[i]]->childrenSize++;
        }

        for(int i = 0; i < n; ++i)
            pArray[i]->children = new NodeR2*[pArray[i]->childrenSize]();

        for(int i = 0; i < n; ++i){
            if(p[i] != -1)
                pArray[p[i]]->children[ptrCounts[p[i]]++] = pArray[i];
        }

        NodeR2 *pR2Root = nullptr;
        for(auto i = 0; i < n; ++i)
            if(p[i] == -1) {
                pR2Root = pArray[i];
                break;
            }

        std::cout << prettyPrint(pR2Root, 0) << "\n";
    }

    if (p != NULL) {
        free(p);
        p = NULL;
    }
    if (repr != NULL) {
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph *graph, Node *start, Node *end, Node *path[]) {
    // TODO: compute the shortest path between the nodes start and end in the given graph
    // the nodes from the path, should be filled, in order, in the array path
    // the number of nodes filled in the path array should be returned
    // if end is not reachable from start, return -1
    // note: the size of the array path is guaranteed to be at least 1000
    if (start == end)
        return 0;

    for (int i = 0; i < graph->nrNodes; ++i)
        graph->v[i]->parent = nullptr;

    bfs(graph, start, NULL);

    int iCount = 0;
    if (end->parent != nullptr) {
        Node* ptrCursor = end;
        while (ptrCursor->parent != nullptr) {
            path[iCount++] = ptrCursor;
            ptrCursor = ptrCursor->parent;
        }
        for (int i = 0; i < iCount / 2; ++i)
            std::swap(path[i], path[iCount - i - 1]);
        return iCount;
    }

    return -1;
}

const int MAX_SIZE = 4500;

void performance() {
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for (n = 1000; n <= MAX_SIZE; n += 100) {
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node **) malloc(graph.nrNodes * sizeof(Node *));
        for (i = 0; i < graph.nrNodes; ++i) {
            graph.v[i] = (Node *) malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }

        for (int i = 0; i < graph.nrNodes; ++i) {
            graph.v[i]->adjSize = 0;
            graph.v[i]->adj = (Node**)malloc(graph.nrNodes * sizeof(Node*));
            memset(graph.v[i]->adj, 0, sizeof(Node*));
        }

        int* ptrEdges = new int[graph.nrNodes]();
        FillRandomArray(ptrEdges, graph.nrNodes, 0, graph.nrNodes - 1, true);
        for (int i = 0; i < graph.nrNodes - 1; ++i)
            graph.v[ptrEdges[i]]->adj[graph.v[ptrEdges[i]]->adjSize++] = graph.v[ptrEdges[i + 1]];


        int k = graph.nrNodes - 1;
        while (k < n) {
            int iV1 = rand() % graph.nrNodes;
            int iV2 = rand() % graph.nrNodes;

            bool bOk = true;

            for (int i = 0; i < graph.nrNodes; ++i)
                for (int j = 0; j < graph.v[i]->adjSize; ++j)
                    if (graph.v[i] == graph.v[iV1] && graph.v[i]->adj[j] == graph.v[iV2] ||
                        graph.v[i] == graph.v[iV2] && graph.v[i]->adj[j] == graph.v[iV1] ||
                        iV1 == iV2)
                        bOk = false;
            if (bOk == true) {
                graph.v[iV1]->adj[graph.v[iV1]->adjSize++] = graph.v[iV2];
                ++k;
            }
        }

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for (n = 100; n <= 200; n += 10) {
        Operation op = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node **) malloc(graph.nrNodes * sizeof(Node *));
        for (i = 0; i < graph.nrNodes; ++i) {
            graph.v[i] = (Node *) malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }

        for (int i = 0; i < graph.nrNodes; ++i) {
            graph.v[i]->adjSize = 0;
            graph.v[i]->adj = (Node**)malloc(n * sizeof(Node*));
            memset(graph.v[i]->adj, 0, sizeof(Node*));
        }

        int* ptrEdges = new int[graph.nrNodes]();
        FillRandomArray(ptrEdges, graph.nrNodes, 0, graph.nrNodes - 1, true);
        for (int i = 0; i < graph.nrNodes - 1; ++i)
            graph.v[ptrEdges[i]]->adj[graph.v[ptrEdges[i]]->adjSize++] = graph.v[ptrEdges[i + 1]];

        int k = graph.nrNodes - 1;
        while (k < 4500) {
            int iV1 = rand() % graph.nrNodes;
            int iV2 = rand() % graph.nrNodes;

            bool bOk = true;

            for (int i = 0; i < graph.nrNodes; ++i)
                for (int j = 0; j < graph.v[i]->adjSize; ++j)
                    if (graph.v[i] == graph.v[iV1] && graph.v[i]->adj[j] == graph.v[iV2] ||
                        graph.v[i] == graph.v[iV2] && graph.v[i]->adj[j] == graph.v[iV1] ||
                        iV1 == iV2)
                        bOk = false;
            if (bOk == true) {
                graph.v[iV1]->adj[graph.v[iV1]->adjSize++] = graph.v[iV2];
                ++k;
            }
        }

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}

Point::Point(int row, int col) {
    this->row = row;
    this->col = col;
}

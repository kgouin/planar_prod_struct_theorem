#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"lca.h"

void LCA_init(struct rmq_struct* s, int** adj, int n){
    //special case
    if (n == 1) return;

    //definitions
    s->n = n;
    s->d = (int*)malloc(((2*(s->n))-1)*sizeof(int));
    s->e = (int*)malloc(((2*(s->n))-1)*sizeof(int));
    s->r = (int*)malloc((s->n)*sizeof(int));
    
    //initialization of representative array
    for (int w = 0; w < (s->n); w++){
        s->r[w] = -1;
    }
    
    //construction of euler tour array, depth array, and representative array
    int u = 0;
    int p = -1;
    int k = 1;
    s->d[0] = 0;
    while(u > -1){
        s->e[k-1] = u;
        if (s->r[u] == -1) s->r[u] = (k-1);
        if (p == adj[u][0]){ //we came from the parent
            p = u;
            if (adj[u][1] > -1){ //left child exists, go there
                u = adj[u][1];
                if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]+1;
            }
            else if (adj[u][2] > -1){ //left child doesn't exist, go to right child
                u = adj[u][2];
                if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]+1;
            }
            else{ //right child doesn't exist, go back to parent
                u = adj[u][0];
                if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]-1;
            }
        }
        else if (p == adj[u][1]){ //we came from left child
            p = u;
            if (adj[u][2] > -1){ //right child exists, go there
                u = adj[u][2];
                if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]+1;
            }
            else { //right child doesn't exist, go back to parent
                u = adj[u][0];
                if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]-1;
            }
        }
        else{ //we came from right child
            //go back to parent
            p = u;
            u = adj[u][0];
            if (k < ((2*(s->n))-1)) s->d[k] = s->d[k-1]-1;
        }
        k++;
    }
    
    RMQ_init(s);
}

int LCA_query(struct rmq_struct* s, int i, int j){
    if (s->r[i] > s->r[j]) return s->e[RMQ_query(s, s->r[j], s->r[i])];
    return s->e[RMQ_query(s, s->r[i], s->r[j])];
}

void LCA_free(struct rmq_struct* s){
    if ((s->n) == 1) return;
    free(s->e);
    free(s->r);
    RMQ_free(s);
}

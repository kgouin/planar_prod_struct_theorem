#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"rmq.h"

int RMQ_ST(struct rmq_struct* s, int i, int j){ // < O(n log n), O(1) >
    //boundary check
    if (i < 0 || j >= (s->n) || j < i) return -1;

    //special case
    if ((s->n) == 1) return 0;

    //sparse table declaration & initialization
    s->st = (int**)calloc((s->n),sizeof(int*));
    for (int i = 0; i < (s->n); i++){
        s->st[i] = (int*)calloc(ceil(log2(s->n)),sizeof(int));
    }

    //sparse table completion
    for (int i = 0; i < (s->n); i++){
        s->st[i][0] = i;
    }
    for (int j = 1; j < (ceil(log2(s->n))); j++){
        for (int i = 0; (i+(1<<j)) <= (s->n); i++){
            if (s->d[s->st[i][j-1]] <= s->d[s->st[i+(1<<(j-1))][j-1]]) s->st[i][j] = s->st[i][j-1];
            else s->st[i][j] = s->st[i+(1<<(j-1))][j-1];
        }
    }

    //perform range-minimum query
    int k, ret;
    ((j-i) == 0) ? (k = 0) : (k = floor(log2(j-i)));
    if (s->d[s->st[i][k]] <= s->d[s->st[j-(1<<k)+1][k]]) ret = s->st[i][k];
    else ret = s->st[j-(1<<k)+1][k];

    return ret;
}

void RMQ_ST_free(struct rmq_struct* s){
    free(s->d);
    if ((s->n) == 1) return;
    for (int i = 0; i < (s->n); i++){
        free(s->st[i]);
    }
    free(s->st);
}

int RMQ_simple(struct rmq_struct* s, int i , int j){
    if (i < 0 || j >= (s->n) || j < i) return -1; //boundary check
    int min = i;
    for (int k = i; k <= j; k++){
        if (s->d[k] < s->d[min]) min = k;
    }
    return min;
}

void RMQ_init(struct rmq_struct* s){
    //special case
    if ((s->n) == 1) return;

    //definitions
    s->b = 2*((8*sizeof(s->n))-__builtin_clz((s->n)-1))/3;
    //define an array A' of size (((s->n)+(s->b)-1)/(s->b)), where A'[i] is the minimum element in the ith block of A
    //define an array B of size (((s->n)+(s->b)-1)/(s->b)), where B[i] is a position in the ith block in which value A'[i] occurs
    s->min_array = (int*)malloc((((s->n)+(s->b)-1)/(s->b)) * sizeof(int));
    s->index_array = (int*)malloc((((s->n)+(s->b)-1)/(s->b)) * sizeof(int));
    int min, min_index;
    for (int k = 0; k < (s->n); k += (s->b)){
        min = s->d[k];
        min_index = 0;
        for (int l = k; (l < (k+(s->b))) && (l < (s->n)); l++){
            if (s->d[l] < min){
                min = s->d[l];
                min_index = l-k;
            }
        }
        s->min_array[k/s->b] = min;
        s->index_array[k/s->b] = min_index;
    }

    //preprocess A' for RMQ
    //sparse table declaration & initialization
    s->st = (int**)calloc((((s->n)+(s->b)-1)/(s->b)),sizeof(int*));
    for (int i = 0; i < (((s->n)+(s->b)-1)/(s->b)); i++){
        s->st[i] = (int*)calloc(ceil(log2(((s->n)+(s->b)-1)/(s->b))),sizeof(int));
    }
    //sparse table completion
    for (int i = 0; i < (((s->n)+(s->b)-1)/(s->b)); i++){
        s->st[i][0] = i;
    }
    for (int j = 1; j < ceil(log2(((s->n)+(s->b)-1)/(s->b))); j++){
        for (int i = 0; (i+(1<<j)) <= (((s->n)+(s->b)-1)/(s->b)); i++){
            if (s->min_array[s->st[i][j-1]] <= s->min_array[s->st[i+(1<<(j-1))][j-1]]) s->st[i][j] = s->st[i][j-1];
            else s->st[i][j] = s->st[i+(1<<(j-1))][j-1];
        }
    }

    //table declaration & initialization
    //t is a master table storing O(√n) tables, each of which has O(b^2) space
    s->t = (int***)calloc(1<<((s->b)-1),sizeof(int**));
    for (int u = 0; u < (1<<((s->b)-1)); u++){
        s->t[u] = (int**)calloc((s->b),sizeof(int*));
        for (int v = 0; v < (s->b); v++){
            s->t[u][v] = (int*)calloc((s->b),sizeof(int));
        }
    }

    //construction of O(b^2) tables
    long z;
    s->signatures = (long*)malloc((((s->n)+(s->b)-1)/(s->b)) * sizeof(long));
    for (int j = 0; j < (((s->n)+(s->b)-1)/(s->b)); j++){
        z = 0;
        for (int l = 0; l < ((s->b)-1); l++){
            if (((((s->b)*j)+l) < (s->n)) && ((((s->b)*j)+l+1) < (s->n)) && (s->d[((s->b)*j)+l] < s->d[((s->b)*j)+l+1])) z += (1<<((s->b)-l-2));
        }
        s->signatures[j] = z;
        if ((s->b) == 1 || s->t[z][1][1] != 1){
            for (int k = 1; k < (s->b); k++){
                s->t[z][k][k] = k;
            }
            for (int k = 0; k < (s->b); k++){
                for (int w = (k+1); w < (s->b); w++){
                    if (((((s->b)*j)+(s->t[z][k][w-1])) <= (s->n)) && ((((s->b)*j)+w) < (s->n)) && (s->d[((s->b)*j)+(s->t[z][k][w-1])] <= s->d[((s->b)*j)+w])) s->t[z][k][w] = s->t[z][k][w-1];
                    else s->t[z][k][w] = w;
                }
            }
        }
    }
}

int RMQ_query(struct rmq_struct* s, int i, int j){
    //boundary check
    if (i < 0 || j >= (s->n) || j < i) return -1;

    //special case
    if ((s->n) == 1) return 0;

    //if i and j are in the same block
    if ((j/(s->b)) == (i/(s->b))) return ((s->b)*(j/(s->b))) + s->t[s->signatures[(j/(s->b))]][i-((j/(s->b))*(s->b))][j-((j/(s->b))*(s->b))];

    //if i and j are in different blocks

    //find min of i to the end of its block
    int suffix_min = ((s->b)*(i/(s->b))) + s->t[s->signatures[(i/(s->b))]][i-(i/(s->b))*(s->b)][(s->b)-1];

    //find the min from j to the beginning of its block
    int prefix_min = ((s->b)*(j/(s->b))) + s->t[s->signatures[(j/(s->b))]][0][j-(j/(s->b))*(s->b)];

    if (((j/(s->b))-1)-((i/(s->b))+1) < 0) return (s->d[suffix_min] <= s->d[prefix_min]) ? suffix_min : prefix_min;

    //find the min of all the blocks in between i's block and j's block
    int range_min, k;
    (((j/(s->b))-1)-((i/(s->b))+1) == 0) ? (k = 0) : (k = floor(log2(((j/(s->b))-1)-((i/(s->b))+1))));
    if (s->min_array[s->st[(i/(s->b))+1][k]] <= s->min_array[s->st[(j/(s->b))-1-(1<<k)+1][k]]) range_min = ((s->st[(i/(s->b))+1][k])*(s->b))+(s->index_array[s->st[(i/(s->b))+1][k]]);
    else range_min = ((s->st[(j/(s->b))-1-(1<<k)+1][k])*(s->b))+(s->index_array[s->st[(j/(s->b))-1-(1<<k)+1][k]]);

    //find true min
    if (s->d[suffix_min] <= s->d[range_min]) return (s->d[suffix_min] <= s->d[prefix_min]) ? suffix_min : prefix_min;
    else return (s->d[range_min] <= s->d[prefix_min]) ? range_min : prefix_min;
}

void RMQ_free(struct rmq_struct* s){
    free(s->d);
    if ((s->n) == 1) return;
    free(s->min_array);
    free(s->index_array);
    for (int i = 0; i < (((s->n)+(s->b)-1)/(s->b)); i++){
        free(s->st[i]);
    }
    free(s->st);
    for (int u = 0; u < (1<<((s->b)-1)); u++){
        for (int v = 0; v < (s->b); v++){
            free(s->t[u][v]);
        }
        free(s->t[u]);
    }
    free(s->t);
    free(s->signatures);
}

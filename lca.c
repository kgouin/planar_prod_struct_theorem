#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

struct rmq_struct{
    int n; //length of array d
    int* e; //array of size (2*n-1) storing the nodes visited in an euler tour of a tree
    int* d; //array of size (2*n-1) storing the depth of the nodes visited in an euler tour of a tree
    int* r; //array of size n storing the representatives of nodes present in e
    int b; //block length
    int* min_array; //array of size (((s->n)+(s->b)-1)/(s->b)) storing the minimum element of each block
    int* index_array; //array of size (((s->n)+(s->b)-1)/(s->b)) storing the in-block index of the minimum element of each block
    int** st; //sparse table of size (((s->n)+(s->b)-1)/(s->b)) by log2((((s->n)+(s->b)-1)/(s->b)))
    int*** t; //master table storing O(√n) tables, each of which has O(b^2) space
    long* signatures; //array of length (((s->n)+(s->b)-1)/(s->b)) storing the integer representation of the binary sequences of the normalized array corresponding to d
};

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
    
    printf("euler tour array:  ");
    for (int k = 0; k < (s->n); k ++){
        printf("%d ", s->e[k]);
    }
    printf("\n");
    printf("depth array:  ");
    for (int k = 0; k < (s->n); k ++){
        printf("%d ", s->d[k]);
    }
    printf("\n");
    printf("representative array:  ");
    for (int k = 0; k < ((s->n)/2)+1; k ++){
        printf("%d ", s->r[k]);
    }
    printf("\n");
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

int main(){
    int k = 10000000;
    int l = 0;
    struct rmq_struct s1;
    clock_t start;
    double elapsed;

    printf("Creating input %d-element input array...", k);
    fflush(stdout);
    start = clock();
    s1.n = k;
    s1.d = (int*)malloc(k * sizeof(int));
    s1.d[0] = 0;
    for (int i = 1; i < k; i++){
        if ((rand() % 2) == 0) s1.d[i] = (s1.d[i-1]) + 1;
        else s1.d[i] = (s1.d[i-1]) - 1;
    }
    printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);

    printf("Creating RMQ structure...");
    fflush(stdout);
    RMQ_init(&s1);
    printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);
    srand(time(NULL));

    printf("Block size = %d\n", s1.b);

    printf("Performing %d queries...", k);
    fflush(stdout);
    while (l < k) {
        int j = (int)(((double)k/RAND_MAX) * rand());
        int i = (int)(((double)(j)/RAND_MAX) * rand());

        RMQ_query(&s1, i, j);
        l++;
    }
    printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);
    RMQ_free(&s1);

    //////////////// LCA testing ////////////////

    int n = 9;
    //create rmq_struct
    struct rmq_struct s;

    //create adjacency list
    int** adj1 = (int**)calloc(n,sizeof(int*));
    for (int i = 0; i < n; i++){
        adj1[i] = (int*)calloc(3,sizeof(int));
    }

    //fill in adjacency list
    //a negative integer represents a null value
    adj1[0][0] = -1; //parent
    adj1[0][1] = 1; //left child
    adj1[0][2] = 5; //right child
    adj1[1][0] = 0;
    adj1[1][1] = 2;
    adj1[1][2] = 3;
    adj1[2][0] = 1;
    adj1[2][1] = -1;
    adj1[2][2] = -1;
    adj1[3][0] = 1;
    adj1[3][1] = -1;
    adj1[3][2] = 4;
    adj1[4][0] = 3;
    adj1[4][1] = -1;
    adj1[4][2] = -1;
    adj1[5][0] = 0;
    adj1[5][1] = 6;
    adj1[5][2] = 7;
    adj1[6][0] = 5;
    adj1[6][1] = -1;
    adj1[6][2] = -1;
    adj1[7][0] = 5;
    adj1[7][1] = 8;
    adj1[7][2] = -1;
    adj1[8][0] = 7;
    adj1[8][1] = -1;
    adj1[8][2] = -1;
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 3; j++){
            if (adj1[i][j] == -1) printf("%d  ", adj1[i][j]);
            else printf("%d   ", adj1[i][j]);
        }
        printf("\n");
    }

    LCA_init(&s, adj1, ((2*(n))-1));
    printf("LCA = %d\n", LCA_query(&s, 2, 1));
    LCA_free(&s);

    //cleanup
    for (int i = 0; i < n; i++){
        free(adj1[i]);
    }
    free(adj1);

    //create adjacency list
    int** adj2 = (int**)calloc(n,sizeof(int*));
    for (int i = 0; i < n; i++){
        adj2[i] = (int*)calloc(3,sizeof(int));
    }

    //fill in adjacency list
    //a negative integer represents a null value
    adj2[0][0] = -1; //parent
    adj2[0][1] = 8; //left child
    adj2[0][2] = 6; //right child
    adj2[1][0] = 6;
    adj2[1][1] = 3;
    adj2[1][2] = -1;
    adj2[2][0] = 8;
    adj2[2][1] = -1;
    adj2[2][2] = 4;
    adj2[3][0] = 1;
    adj2[3][1] = -1;
    adj2[3][2] = -1;
    adj2[4][0] = 2;
    adj2[4][1] = -1;
    adj2[4][2] = -1;
    adj2[5][0] = 8;
    adj2[5][1] = -1;
    adj2[5][2] = -1;
    adj2[6][0] = 0;
    adj2[6][1] = 7;
    adj2[6][2] = 1;
    adj2[7][0] = 6;
    adj2[7][1] = -1;
    adj2[7][2] = -1;
    adj2[8][0] = 0;
    adj2[8][1] = 5;
    adj2[8][2] = 2;

    printf("----------------\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 3; j++){
            if (adj2[i][j] == -1) printf("%d  ", adj2[i][j]);
            else printf("%d   ", adj2[i][j]);
        }
        printf("\n");
    }

    LCA_init(&s, adj2, ((2*(n))-1));
    printf("LCA = %d\n", LCA_query(&s, 2, 1));
    LCA_free(&s);

    //cleanup
    for (int i = 0; i < n; i++){
        free(adj2[i]);
    }
    free(adj2);

 //create adjacency list
    int** adj3 = (int**)calloc(n,sizeof(int*));
    for (int i = 0; i < n; i++){
        adj3[i] = (int*)calloc(3,sizeof(int));
    }

    //fill in adjacency list
    //a negative integer represents a null value
    adj3[0][0] = -1; //parent
    adj3[0][1] = 5; //left child
    adj3[0][2] = 7; //right child
    adj3[1][0] = 5;
    adj3[1][1] = -1;
    adj3[1][2] = 6;
    adj3[2][0] = 5;
    adj3[2][1] = -1;
    adj3[2][2] = -1;
    adj3[3][0] = 7;
    adj3[3][1] = 8;
    adj3[3][2] = -1;
    adj3[4][0] = 7;
    adj3[4][1] = -1;
    adj3[4][2] = -1;
    adj3[5][0] = 0;
    adj3[5][1] = 2;
    adj3[5][2] = 1;
    adj3[6][0] = 1;
    adj3[6][1] = -1;
    adj3[6][2] = -1;
    adj3[7][0] = 0;
    adj3[7][1] = 4;
    adj3[7][2] = 3;
    adj3[8][0] = 3;
    adj3[8][1] = -1;
    adj3[8][2] = -1;

    printf("----------------\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 3; j++){
            if (adj3[i][j] == -1) printf("%d  ", adj3[i][j]);
            else printf("%d   ", adj3[i][j]);
        }
        printf("\n");
    }

    LCA_init(&s, adj3, ((2*(n))-1));
    printf("LCA = %d\n", LCA_query(&s, 2, 1));
    LCA_free(&s);

    //cleanup
    for (int i = 0; i < n; i++){
        free(adj3[i]);
    }
    free(adj3);
}

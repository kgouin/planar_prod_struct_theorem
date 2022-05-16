#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int RMQ_ST(int* in_array, int i, int j, int n){ // < O(n log n), O(1) >

    //boundary check
    if (i < 0 || j >= n) exit(0);

    //sparse table declaration & initialization
    int** st;
    st = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++){
        st[i] = (int*)malloc(log2(n) * sizeof(int));
    }

    //sparse table completion
    for (int i = 0; i < n; i++){
        st[i][0] = in_array[i];
    }
    for (int j = 1; j <= log2(n); j++){
        for (int i = 0; (i+(1<<j)) <= n; i++){
            if (st[i][j-1] <= st[i+(1<<(j-1))][j-1]) st[i][j] = st[i][j-1];
            else st[i][j] = st[i+(1<<(j-1))][j-1];
        }
    }

    //sparse table print
    for (int i = 0; i < n; i++){
        for (int j = 0; j <= log2(n); j++){
            printf("%d  ", st[i][j]);
        }
        printf("\n");
    }
    
    //perform range-minimum query
    int ret;
    int k = floor(log2(j-i));
    if (st[i][k] <= st[j-(1<<k)+1][k]) ret = st[i][k];
    else ret = st[j-(1<<k)+1][k];

    //cleanup
    for (int i = 0; i < n; i++){
        free(st[i]);
    }
    free(st);

    //find index of minimum value (not sure about this in terms of efficiency)
    int m = i;
    while(m < n){
        if (in_array[m] == ret) return m;
        m++;
    }

    exit(0); //somehow the minimum element we found was not in our array. (we should never reach this)

}


int PlusMinusOne_RMQ(int* in_array, int i, int j, int n){ // < O(n), O(1) >
    //let's assume for now that log2(n) produces an integer, and assume an array full of complete blocks

    //boundary check
    if (i < 0 || j >= n) exit(0);

    //definitions
    int numBlocks = (int)(n/(log2(n)/2));
    int b = (floor(log2(n)))/2;
    float specific_block_i = (i/((floor(log2(n)))/2))+1;
    float specific_block_j = (j/((floor(log2(n)))/2))+1;
    int block_i = (int)((i/((floor(log2(n)))/2))+1); //assuming that the first block is B1 (and not B0)
    int block_j = (int)((j/((floor(log2(n)))/2))+1); //assuming that the first block is B1
    int block_i_end = (block_i*b)-1;
    int block_j_start = (block_j*b)-b;

    //testing
    /*printf("n = %d\n", n);
    printf("i = %d\n", i);
    printf("j = %d\n", j);
    printf("b = %d\n", b);
    printf("there are %d blocks\n", numBlocks);
    printf("specific block i = %f\n", specific_block_i);
    printf("specific block j = %f\n", specific_block_j);
    printf("i is in block B%d\n", block_i);
    printf("j is in block B%d\n", block_j);
    printf("the block which contains i ends at index %d\n", block_i_end);
    printf("the block which contains j starts at index %d\n", block_j_start);*/

    //testing
    printf("in_array before normalization\n");
    for (int m = 0; m < n; m++){
        printf("%d  ", in_array[m]);
    }
    printf("\n");

    //normalize all blocks
    int initialOffset;
    for (int k = 0; k < n; k += b){
        initialOffset = in_array[k];
        for (int l = k; l < k+b; l++){
            in_array[l] -= initialOffset;
        }
    }

    //testing
    printf("in_array after normalization\n");
    for (int m = 0; m < n; m++){
        printf("%d  ", in_array[m]);
    }
    printf("\n");

    //create 0(√n) tables, one for each possible normalized block
    //in each table, we put all b^2 answers to all in-block queries (this is what we will be using below)

    //if i and j are in different blocks:
    /*if (block_i != block_j){
        //find min of i to the end of its block
        int suffixMin = in_array[i];
        for (int k = i; k < block_i_end; k++){ //temporary. not the correct approach (see above)
            if (in_array[k] < suffixMin) suffixMin = in_array[k]; //temporary. not the correct approach (see above)
        }
        //find the min of all the blocks in between i's block and j's block
        //...
        //find the min from j to the beginning of its block
        int prefixMin = in_array[j];
        for (int k = j; k > block_j_start; k--){ //temporary. not the correct approach (see above)
            if (in_array[k] < suffixMin) suffixMin = in_array[k]; //temporary. not the correct approach (see above)
        }
    }*/

    return 0;

}


int main(){

    int test1[] = {0,1,2,1,3,2,1,2,3,2,3,4,3,2,1,0};
    int test2[] = {0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,-1};
    int min1 = PlusMinusOne_RMQ(test1, 2, 9, (sizeof(test1)/sizeof(*test1)));
    //int min2 = PlusMinusOne_RMQ(test2, 13, 254, (sizeof(test2)/sizeof(*test2)));
    //printf("min1 is at index %d\nmin2 is at index %d\n", min1, min2);

}
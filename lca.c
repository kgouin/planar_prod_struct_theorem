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

    //find index of minimum value (not the correct approach)
    //we will probably need to use an additional array of size n to keep track of where the minima came from
    int m = i;
    while(m < n){
        if (in_array[m] == ret) return m;
        m++;
    }

    exit(0); //somehow the minimum element we found was not in our array (we should never reach this)

}


int PlusMinusOne_RMQ(int* in_array, int i, int j, int n){ // < O(n), O(1) >
    //let's assume for now that log2(n) produces an integer, and assume an array full of complete blocks

    //boundary check
    if (i < 0 || j >= n) exit(0);

    //definitions (assuming that the first block is B0)
    int num_blocks = (int)(n/(log2(n)/2));
    int b = (floor(log2(n)))/2;
    float specific_block_i = (i/((floor(log2(n)))/2));
    float specific_block_j = (j/((floor(log2(n)))/2));
    int block_i = (int)((i/((floor(log2(n)))/2)));
    int block_j = (int)((j/((floor(log2(n)))/2)));
    int block_i_end = ((block_i+1)*b)-1;
    int block_j_start = ((block_j+1)*b)-b;

    //testing
    /*printf("n = %d\n", n);
    printf("i = %d\n", i);
    printf("j = %d\n", j);
    printf("b = %d\n", b);
    printf("there are %d blocks\n", num_blocks);
    printf("specific block i = %f\n", specific_block_i);
    printf("specific block j = %f\n", specific_block_j);
    printf("i is in block B%d\n", block_i);
    printf("j is in block B%d\n", block_j);
    printf("the block which contains i ends at index %d\n", block_i_end);
    printf("the block which contains j starts at index %d\n", block_j_start);*/

    //define an array A' of size num_blocks, where A'[i] is the minimum element in the ith block of A
    //define an array B of size num_blocks, where B[i] is a position in the ith block in which value A'[i] occurs
    int* min_array = (int*)malloc(num_blocks * sizeof(int));
    int* index_array = (int*)malloc(num_blocks * sizeof(int));
    int min, min_index;
    for (int k = 0; k < n; k += b){ // O(n) time
        min = in_array[k];
        min_index = 0;
        for (int l = k; l < k+b; l++){
            if (in_array[l] < min){
                min = in_array[l];
                min_index = l-k;
            }
        }
        min_array[k/b] = min;
        index_array[k/b] = min_index;
    }

    //testing
    printf("min_array:    ");
    for (int k = 0; k < num_blocks; k++){
        printf("%d  ", min_array[k]);
    }
    printf("\n");
    printf("index_array:  ");
    for (int k = 0; k < num_blocks; k++){
        printf("%d  ", index_array[k]);
    }
    printf("\n");

    //preprocess A' for RMQ
    //sparse table declaration & initialization
    int** st;
    st = (int**)malloc(num_blocks * sizeof(int*));
    for (int i = 0; i < num_blocks; i++){
        st[i] = (int*)malloc(log2(num_blocks) * sizeof(int));
    }
    //sparse table completion
    for (int i = 0; i < num_blocks; i++){
        st[i][0] = min_array[i];
    }
    for (int j = 1; j <= log2(num_blocks); j++){
        for (int i = 0; (i+(1<<j)) <= num_blocks; i++){
            if (st[i][j-1] <= st[i+(1<<(j-1))][j-1]) st[i][j] = st[i][j-1];
            else st[i][j] = st[i+(1<<(j-1))][j-1];
        }
    }
    //sparse table print
    for (int i = 0; i < num_blocks; i++){
        for (int j = 0; j <= log2(num_blocks); j++){
            printf("%d  ", st[i][j]);
        }
        printf("\n");
    }

    //normalize all blocks (this might have to be a binary array)
    int initialOffset;
    for (int k = 0; k < n; k += b){ // O(n) time
        initialOffset = in_array[k];
        for (int l = k; l < k+b; l++){
            in_array[l] -= initialOffset;
        }
    }

    //testing
    printf("normalized in_array: ");
    for (int m = 0; m < n; m++){
        printf("%d  ", in_array[m]);
    }
    printf("\n");

    //for in-block queries:
    //create O(√n) tables, one for each possible kind of normalized block
    //in each table, we put all b^2 answers to all in-block queries (this is what we will be using below)
    //...

    //if i and j are in the same block:
    //use one of O(√n) normalized block tables
    //...

    //if i and j are in different blocks:
    if (block_i != block_j){
        //find min of i to the end of its block (O(1) time)
        //use one of O(√n) normalized block tables
        //...

        //find the min of all the blocks in between i's block and j's block (O(1) time)
        //this works. do not modify
        int range_min;
        int k = floor(log2((block_j-1)-(block_i+1)));
        if (st[block_i+1][k] <= st[block_j-1-(1<<k)+1][k]) range_min = st[block_i+1][k];
        else range_min = st[block_j-1-(1<<k)+1][k];
        printf("range_min = %d\n", range_min);

        //find the min from j to the beginning of its block (O(1) time)
        //use one of O(√n) normalized block tables
        //...

    }

    //cleanup
    for (int i = 0; i < num_blocks; i++){
        free(st[i]);
    }
    free(st);
    free(min_array);
    free(index_array);

    return 0;

}


int main(){

    int test1[] = {0,1,2,1,3,2,1,2,3,2,3,4,3,2,1,0};
    int test2[] = {0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,-1};
    int min1 = PlusMinusOne_RMQ(test1, 2, 9, (sizeof(test1)/sizeof(*test1)));
    //int min2 = PlusMinusOne_RMQ(test2, 13, 254, (sizeof(test2)/sizeof(*test2)));
    //printf("min1 is at index %d\nmin2 is at index %d\n", min1, min2);

}
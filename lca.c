#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//struct for RMQ

//RMQ init
//RMQ query
//RMQ free


int RMQ_ST(int* in_array, int i, int j, int n){ // < O(n log n), O(1) >

    //boundary check
    if (i < 0 || j >= n || j < i) return -1;

    //sparse table declaration & initialization
    int** st;
    st = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++){
        st[i] = (int*)malloc(log2(n) * sizeof(int));
    }

    //sparse table completion
    for (int i = 0; i < n; i++){
        st[i][0] = i;
    }
    for (int j = 1; j <= log2(n); j++){
        for (int i = 0; (i+(1<<j)) <= n; i++){
            if (in_array[st[i][j-1]] <= in_array[st[i+(1<<(j-1))][j-1]]) st[i][j] = st[i][j-1];
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
    if (in_array[st[i][k]] <= in_array[st[j-(1<<k)+1][k]]) ret = st[i][k];
    else ret = st[j-(1<<k)+1][k];

    //cleanup
    for (int i = 0; i < n; i++){
        free(st[i]);
    }
    free(st);

    return ret;

}


int PlusMinusOne_RMQ(int* in_array, int i, int j, int n){ // < O(n), O(1) >

    //boundary check
    if (i < 0 || j >= n || j < i) return -1;

    //definitions (assuming that the first block is B0)
    int b = ceil((double)(log2(n)/2));
    int num_blocks = ceil(n/(double)b);
    int block_i = (i/b);
    int block_j = (j/b);
    int block_i_end = ((block_i+1)*b)-1;
    int block_j_start = ((block_j+1)*b)-b;

    //testing
    /*printf("n = %d\n", n);
    printf("i = %d\n", i);
    printf("j = %d\n", j);
    printf("b = %d\n", b);
    printf("there are %d blocks\n", num_blocks);
    printf("i is in block B%d\n", block_i);
    printf("j is in block B%d\n", block_j);
    printf("the block which contains i ends at index %d\n", block_i_end);
    printf("the block which contains j starts at index %d\n", block_j_start);*/

    //define an array A' of size num_blocks, where A'[i] is the minimum element in the ith block of A
    //define an array B of size num_blocks, where B[i] is a position in the ith block in which value A'[i] occurs
    //do we even need array A' ?
    int* min_array = (int*)malloc(num_blocks * sizeof(int)); //our array A'
    int* index_array = (int*)malloc(num_blocks * sizeof(int)); //our array B
    int min, min_index;
    for (int k = 0; k < n; k += b){
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

    //preprocess A' for RMQ
    //sparse table declaration & initialization
    int** st;
    st = (int**)malloc(num_blocks * sizeof(int*));
    for (int i = 0; i < num_blocks; i++){
        st[i] = (int*)malloc(log2(num_blocks) * sizeof(int));
    }
    //sparse table completion
    for (int i = 0; i < num_blocks; i++){
        st[i][0] = i;
    }
    for (int j = 1; j <= log2(num_blocks); j++){
        for (int i = 0; (i+(1<<j)) <= num_blocks; i++){
            if (min_array[st[i][j-1]] <= min_array[st[i+(1<<(j-1))][j-1]]) st[i][j] = st[i][j-1];
            else st[i][j] = st[i+(1<<(j-1))][j-1];
        }
    }

    //table declaration & initialization
    //t is a master table containing O(√n) tables, each of which has O(b^2) space
    int *** t = (int***)malloc(sqrt(n) * sizeof(int**));
    for (int u = 0; u < sqrt(n); u++){
        t[u] = (int**)malloc((b+1) * sizeof(int*));
        for (int v = 0; v < (b+1); v++){
            t[u][v] = (int*)malloc((b+1) * sizeof(int));
        }
    }

    //construction of O(b^2) tables
    char* temp; //maybe do the string to long converion myself so that I don't need this temp variable
    long z;
    char** normalized_array = (char**)malloc((num_blocks) * sizeof(char));
    for (int v = 0; v < num_blocks; v++){
        normalized_array[v] = (char*)malloc((b-1) * sizeof(char));
    }
    long* signatures = (long*)malloc(num_blocks * sizeof(long));
    for (int j = 0; j < num_blocks; j++){
        for (int l = 0; l < b-1; l++){
            if (((b*j)+l) < n && ((b*j)+l+1) < n && in_array[(b*j)+l] < in_array[(b*j)+l+1]) normalized_array[j][l] = '1';
            else normalized_array[j][l] = '0';
        }
        z = strtol(normalized_array[j], &temp, 2); //store integer representation of binary sequence in signatures array
        signatures[j] = z;
        for (int k = 1; k < b+1; k++){
            t[z][k][k] = k;
        }
        for (int k = 0; k < b+1; k++){
            for (int w = k+1; w < b+1; w++){
                if ((b*j)+(t[z][k][w-1]) < n && ((b*j)+w) < n && in_array[(b*j)+(t[z][k][w-1])] < in_array[(b*j)+w]) t[z][k][w] = t[z][k][w-1];
                else t[z][k][w] = w;
            }
        }
    }

    //if i and j are in the same block
    if (block_i == block_j){
        int short_min = (b*block_j) + t[signatures[block_j]][i-block_j_start][j-block_j_start];
        printf("short_min = %d\n", short_min);
        return short_min;
    }

    //if i and j are in different blocks
    if (block_i != block_j){
        //find min of i to the end of its block
        int suffix_min = (b*block_i) + t[signatures[block_i]][i-(block_i_end-b+1)][b-1];
        printf("suffix_min = %d\n", suffix_min);

        //find the min of all the blocks in between i's block and j's block
        int range_min, k;
        ((block_j-1)-(block_i+1) == 0) ? (k = 0) : (k = floor(log2((block_j-1)-(block_i+1))));

        if (min_array[st[block_i+1][k]] <= min_array[st[block_j-1-(1<<k)+1][k]]){
            //range_min = index_array[st[block_i+1][k]];
            //printf("range_min is in block %d, index %d\n", st[block_i+1][k], range_min);
            //printf("range_min is in block starting at index %d\n", ((st[block_i+1][k])*b));
            //printf("range_min = %d\n", ((st[block_i+1][k])*b)+(index_array[st[block_i+1][k]]));
            range_min = ((st[block_i+1][k])*b)+(index_array[st[block_i+1][k]]);
        }
        else{
            //range_min = index_array[st[block_j-1-(1<<k)+1][k]];
            //printf("range_min is in block %d, index %d\n", st[block_j-1-(1<<k)+1][k], range_min);
            //printf("range_min is in block starting at index %d\n", ((st[block_j-1-(1<<k)+1][k])*b));
            //printf("range_min = %d\n", ((st[block_j-1-(1<<k)+1][k])*b)+(index_array[st[block_j-1-(1<<k)+1][k]]));
            range_min = ((st[block_j-1-(1<<k)+1][k])*b)+(index_array[st[block_j-1-(1<<k)+1][k]]);
        }

        printf("range_min = %d\n", range_min);

        //find the min from j to the beginning of its block
        int prefix_min = (b*block_j) + t[signatures[block_j]][0][j-block_j_start];
        printf("prefix_min = %d\n", prefix_min);

        //find ultimate min
        if (in_array[suffix_min] <= in_array[range_min]) return (in_array[suffix_min] <= in_array[prefix_min]) ? suffix_min : prefix_min;
        else return (in_array[range_min] <= in_array[prefix_min]) ? range_min : prefix_min;

    }

    //cleanup
    /*for (int i = 0; i < num_blocks; i++){
        free(st[i]);
    }
    free(st);
    for (int u = 0; u < sqrt(n); u++){
        for (int v = 0; v < b; v++){
            free(t[u][v]);
        }
        free(t[u]);
    }
    free(t);
    free(min_array);
    free(index_array);
    free(signatures);*/

    return -1; //error

}


int main(){

    int test1[] = {-1,0,1,2,3,2,1,2,3,2,3,4,3,2,1,0};
    int test2[] = {0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,-1};
    int test3[] = {0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,10,9,8,7,6,7,8,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,0,1,2,3,4,5,6,5,4,3,2,1,0,-1,-2};
    //int min = RMQ_ST(test1, 4, 8, sizeof(test1)/sizeof(*test1));
    //int min = RMQ_ST(test2, 2, 42, sizeof(test2)/sizeof(*test2));
    //int min = PlusMinusOne_RMQ(test1, 1, 9, (sizeof(test1)/sizeof(*test1)));
    //int min = PlusMinusOne_RMQ(test2, 18, 42, (sizeof(test2)/sizeof(*test2)));
    //int min = PlusMinusOne_RMQ(test2, 1, 9, (sizeof(test2)/sizeof(*test2)));
    //int min = PlusMinusOne_RMQ(test2, 2, 3, (sizeof(test1)/sizeof(*test1)));
    //int min = PlusMinusOne_RMQ(test2, 16, 18, (sizeof(test2)/sizeof(*test2)));
    //int min = PlusMinusOne_RMQ(test3, 18, 42, (sizeof(test3)/sizeof(*test3)));
    //int min = PlusMinusOne_RMQ(test3, 2, 254, (sizeof(test3)/sizeof(*test3)));
    int min = PlusMinusOne_RMQ(test3, 18, 254, (sizeof(test3)/sizeof(*test3)));
    printf("min is at index %d\n", min);

}
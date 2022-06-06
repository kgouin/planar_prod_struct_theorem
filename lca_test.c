#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"lca.h"

int main(){
    int n = 9;
    //create lca_struct
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

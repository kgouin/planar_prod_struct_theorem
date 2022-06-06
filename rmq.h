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

void RMQ_init(struct rmq_struct*);
int RMQ_query(struct rmq_struct*, int, int);
void RMQ_free(struct rmq_struct*);

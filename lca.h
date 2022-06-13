#include"rmq.h"

int LCA_simple(int**, int, int, int);
void LCA_init(struct rmq_struct*, int**, int);
int LCA_query(struct rmq_struct*, int, int);
void LCA_free(struct rmq_struct*);

#include"lca.h"
#include"bfs.h"

struct tripod_decomposition_struct{
	int v_x[3];
	int v_x_next[3];
	int v_x_op[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	int* vertex_tripod_assign; //for each vertex v, vertex_tripod_assign[v] represents the colour of the tripod v belongs to, -1 if it does not yet belong to a tripod
	int* face_tripod_assign; //for each face f, face_tripod_assign[f] is -1 if it does not belong to a tripod, and face_tripod_assign[f]=f otherwise
	int* tripod_assign_order; //for each tripod, tripod_assign_order stores its sp colour in the order in which we find them
	int tripod_assign_order_index; //index at which we add tripod to t->tripod_assign_order
};

void init(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*);
int* decompose(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*);

void store_tripod(struct bfs_struct*, struct tripod_decomposition_struct*, int);

int* trichromatic_tripod(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int, int, int);
void trichromatic_orient_subproblems(struct bfs_struct*, struct tripod_decomposition_struct*, int, int, int, int);
void trichromatic_decompose(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int, int, int, int);

int* bichromatic_tripod(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int, int);
void bichromatic_decompose(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int, int, int);

int* monochromatic_tripod(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int);
void monochromatic_decompose(struct bfs_struct*, struct rmq_struct*, struct tripod_decomposition_struct*, int, int);

void tprint(struct tripod_decomposition_struct*);

void tripod_free(struct tripod_decomposition_struct*);

int three_tree_test(struct bfs_struct*, struct tripod_decomposition_struct*, int);
int three_tree_test_pt2(struct bfs_struct* b, struct tripod_decomposition_struct* t);
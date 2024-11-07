#include"lca.h"
#include"bfs.h"

struct tripod_decomposition_struct{
	int v_a; //the last vertex before we hit another tripod
	int v_a_next; //the vertex at which we hit another tripod
	int v_a_op; //the face with bichromatic edge on the cycle defining the subproblem
	int v_a_mirror; //the face adjacent to newly found sp
	int v_b;
	int v_b_next;
	int v_b_op;
	int v_b_mirror;
	int v_c;
	int v_c_next;
	int v_c_op;
	int v_c_mirror;
	int* v_x;
	int* v_x_next;
	int* v_x_op;
	int* v_x_mirror;
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
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

int main(int argc, char **argv){
	if (argc < 5){
		printf("Usage:\n");
		printf("This program requires 5 arguments, organized as follows:\n");
		printf("argv[0] is the name of the executable\n");
		printf("argv[1] is the name of the adjacencies text file\n");
		printf("argv[2] is the name of the simplicies text file\n");
		printf("argv[3] is the name of the triangles text file\n");
		printf("argv[4] is the basename for output text files\n");
		return 0;
	}

	struct bfs_struct b;
	struct rmq_struct r;
	struct tripod_decomposition_struct t;

	init(&b, &r, &t, argv[1], argv[2], argv[3]);
	decompose(&b, &r, &t);

	for (int k = 0; k < ((&b)->v); k++){
		if (((&t)->vertex_tripod_assign)[k] == -1){
			printf("error. some vertices are not labeled as belonging to a tripod.\n");
			exit(0);
		}
	}

	for (int k = 0; k < ((&b)->f); k++){
		if (((&t)->face_tripod_assign)[k] == -1){
			printf("error. some faces are not labeled as belonging to a tripod.\n");
			exit(0);
		}
	}

	if (!three_tree_test(&b, &t)) printf("three tree test unsuccessful\n");

	//write bfs tree to file (1 int per line, b.v lines)
	char* complete_bfs_tree_file_name;
	asprintf(&complete_bfs_tree_file_name, "%s_bfs_tree.txt", argv[4]);
	FILE *fd;
	fd = fopen(complete_bfs_tree_file_name, "w");
	if (!fd) exit(0);
	fprintf(fd, "%d\n", b.v);
	for (int i = 0; i < b.v; i++){
		fprintf(fd, "%d\n", b.bt[i]);
	}
	fclose(fd);

	//write vertex_tripod_assign array to file (1 int per line, b.v lines)
	char* complete_vertex_tripod_assign_file_name;
	asprintf(&complete_vertex_tripod_assign_file_name, "%s_vertex_tripod_assign.txt", argv[4]);
	fd = fopen(complete_vertex_tripod_assign_file_name, "w");
	if (!fd) exit(0);
	fprintf(fd, "%d\n", b.v);
	for (int i = 0; i < b.v; i++){
		fprintf(fd, "%d\n", t.vertex_tripod_assign[i]);
	}
	fclose(fd);

	//write tripod_adjacency_list array to file (3 ints per line, b.f lines)
	char* complete_tripod_adjacency_list_file_name;
	asprintf(&complete_tripod_adjacency_list_file_name, "%s_tripod_adjacency_list.txt", argv[4]);
	fd = fopen(complete_tripod_adjacency_list_file_name, "w");
	if (!fd) exit(0);
	fprintf(fd, "%d\n", b.f);
	for (int i = 0; i < b.f; i++){
		fprintf(fd, "%d %d %d\n", t.tripod_adjacency_list[i][0], t.tripod_adjacency_list[i][1], t.tripod_adjacency_list[i][2]);
	}
	fclose(fd);

	//we could have one more file, where, for each tripod, for each leg of given tripod, we list its vertices following the upward bfs path
	//tripods -> groups of three lines with many ints (vertices) per line, tripod_assign_order_index many lines

	tripod_free(&b, &t);
	LCA_free(&r);
	BFS_free(&b);
}

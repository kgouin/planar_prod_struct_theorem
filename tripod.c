#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void init(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){
	//initialize bfs_struct & rmq_struct
	BFS_init(b);
	BFS(b);
	LCA_init(r, b->ct, (2*(b->f))-1);

	//initialize tripod_decomposition_struct
	t->vertex_tripod_assign = malloc((b->v)*sizeof(int));
	//label the three vertices incident to face 0 as belonging to three exterior tripods
	t->vertex_tripod_assign[b->sim[0][0]] = (b->f);
	t->vertex_tripod_assign[b->sim[0][1]] = (b->f)+1;
	t->vertex_tripod_assign[b->sim[0][2]] = (b->f)+2;
	for (int i = 3; i < (b->v); i++){
		t->vertex_tripod_assign[i] = -1; //label remaining vertices with -1
	}
	t->face_tripod_assign = malloc((b->f)*sizeof(int));
	t->face_tripod_assign[0] = 0; //this is a special case, since 0 is the outer face
	for (int i = 1; i < (b->f); i++){
		t->face_tripod_assign[i] = -1;
	}
	t->tripod_assign_order = malloc(((b->f)+2)*sizeof(int));
	//label the first three entries as the three exterior tripods that constitute our base case
	t->tripod_assign_order[0] = (b->f);
	t->tripod_assign_order[1] = (b->f)+1;
	t->tripod_assign_order[2] = (b->f)+2;
	for (int i = 3; i < ((b->f)+2); i++){
		t->tripod_assign_order[i] = -1;
	}
	t->tripod_assign_order_index = 3;

	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("initialization complete\n");

	printf("bfs tree = "); //testing
	for (int i = 0; i < b->v; i++){
		printf("%d ", b->bt[i]);
	}
	printf("\n");

	//write bfs to file for later visualization
	FILE *fd;
	fd = fopen("bfs.txt", "w");
	if (!fd) exit(1);
	fprintf(fd, "%d\n", b->v);
	for (int i = 0; i < b->v; i++){
		fprintf(fd, "%d\n", b->bt[i]);
	}
	fclose(fd);
}

int* decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){

	//start decomposition with the three triangles adjacent to outer face
	trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1]); //these are in counterclockwise order

	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("decomposition complete\n");

	//t->face_tripod_assign and t->tripod_assign_order exist for testing purposes
	printf("face_tripod_assign = [ ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("]\n");
	printf("tripod_assign_order = [ ");
	for (int i = 0; i < t->tripod_assign_order_index; i++){
		printf("%d ", t->tripod_assign_order[i]);
	}
	printf("]\n");

	return t->vertex_tripod_assign;
}

/******************************************************** trichromatic ********************************************************/

int* trichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3){
	//sperner triangle identification
	int sp;
	if (f1 == f2 && f1 == f3){
		sp = f1;
		t->face_tripod_assign[sp] = sp; //keep track of sperner triangles
		//don't add tripod to t->tripod_assign_order
		return 0;
	}
	else if (f1 != f2 && f1 != f3 && f2 != f3){
		//find sperner triangle
		if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3) && LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		else {
			if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3)) sp = LCA_query(r, f2, f3);
			else if (LCA_query(r, f1, f2) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f3);
			else if (LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		}
	}
	else { //if two of the three input triangles are the same
		//then they are (it is) the sperner triangle
		//performing LCA on f1, f2, and f3 may lead to an incorrect selection of sp
		if (f1 == f2) sp = f1;
		else sp = f3;
	}

	//store tripod
	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_x[0]] == sp || t->vertex_tripod_assign[t->v_x[1]] == sp || t->vertex_tripod_assign[t->v_x[2]] == sp){ //if one of our legs is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[f1] == f1 || t->face_tripod_assign[f2] == f2 || t->face_tripod_assign[f3] == f3){
			printf("one of our input triangles is a sperner triangle. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
			exit(0);
		}
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	if (f1 == sp && f2 == sp && f3 == sp) return 0; //no subproblems exist

	trichromatic_orient_subproblems(b, t, sp, f1, f2, f3);

	tprint(t);

	trichromatic_decompose(b, r, t, sp, f1, f2, f3);

	return 0;
}

void trichromatic_orient_subproblems(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3){

	int leg_colours[3];
	for (int i = 0; i < 3; i++){ //identify the endpoint colour of each leg
		(t->v_x_next[i] == -1) ? (leg_colours[i] = t->vertex_tripod_assign[t->v_x[i]]) : (leg_colours[i] = t->vertex_tripod_assign[t->v_x_next[i]]);
	}

	int f[3] = {f1, f2, f3};

	int k;
	for (k = 0; k < 3 && f[k] == sp; k++){
		//this loop checks which subproblem is non-empty
	}
	//use: f[k], f[(k+1)%3], f[(k+2)%3]
	//rotate f

	printf("f = f%d = %d\n", k+1, f[k]); //temporary, so that diff will show no difference between the original code and the new modular arithmetic code

	int double_match = 0;
	for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
		if (t->vertex_tripod_assign[b->sim[f[k]][i]] == leg_colours[0]){
			for (int j = 0; j < 3; j++){
				if (t->vertex_tripod_assign[b->sim[f[k]][j]] == leg_colours[1]){
					double_match = 1;
				}
			}
		}
	}
	if (double_match){
		printf("ORIENTATION 1\n");
		t->v_x_op[0] = f[k];
		t->v_x_op[1] = f[(k+1)%3];
		t->v_x_op[2] = f[(k+2)%3];
	}
	else { //no double match. reset values and try with leg_colour_b and leg_colour_c
		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
			if (t->vertex_tripod_assign[b->sim[f[k]][i]] == leg_colours[1]){
				for (int j = 0; j < 3; j++){
					if (t->vertex_tripod_assign[b->sim[f[k]][j]] == leg_colours[2]){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("ORIENTATION 2\n");
			t->v_x_op[0] = f[(k+2)%3];
			t->v_x_op[1] = f[k];
			t->v_x_op[2] = f[(k+1)%3];
		}
		else { //no double match. f must match with leg_colour_c and leg_colour_a.
			printf("ORIENTATION 3\n");
			t->v_x_op[0] = f[(k+1)%3];
			t->v_x_op[1] = f[(k+2)%3];
			t->v_x_op[2] = f[k];
		}
	}
}

void trichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3){
	//definitions
	int v_x[3];
	int v_x_next[3];
	int v_x_op[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_op[i] = t->v_x_op[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	for (int i = 0; i < 3; i++){
		if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is non-empty && leg (i+1)%3 is non-empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_r[(i+1)%3], v_x_l[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is non-empty && leg (i+1)%3 is empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_mirror[i], v_x_l[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_mirror[i], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is empty && leg (i+1)%3 is non-empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_r[(i+1)%3], v_x_mirror[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_mirror[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is empty && leg (i+1)%3 is empty
			if (t->face_tripod_assign[v_x_op[i]] > -1){
				printf("no subproblem x for sp %d\n", sp);
				v_x_mirror[i] = -1; //set v_x_mirror[i] to negative value to later check whether it is in our subproblem boundary
			}
			else {
				printf("subproblem x for sp %d is bichromatic\n", sp);
				printf("subproblem x for sp %d will be on faces %d, %d\n", sp, v_x_op[i], v_x_mirror[i]);
				bichromatic_tripod( b, r, t, v_x_op[i], v_x_mirror[i]);
			}
		}
	}
}

/******************************************************** bichromatic ********************************************************/

int* bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2){
	int sp;

	sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_x[0]] == sp || t->vertex_tripod_assign[t->v_x[1]] == sp || t->vertex_tripod_assign[t->v_x[2]] == sp){ //if one of our legs is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[f1] == f1 || t->face_tripod_assign[f2] == f2){
			printf("one of our input triangles is a sperner triangle. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
			exit(0);
		}
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	printf("v_a = %d\n", t->v_x[0]);
	printf("v_b = %d\n", t->v_x[1]);
	printf("v_c = %d\n", t->v_x[2]);
	printf("v_a_next = %d\n", t->v_x_next[0]);
	printf("v_b_next = %d\n", t->v_x_next[1]);
	printf("v_c_next = %d\n", t->v_x_next[2]);
	printf("v_a_mirror = %d\n", t->v_x_mirror[0]);
	printf("v_b_mirror = %d\n", t->v_x_mirror[1]);
	printf("v_c_mirror = %d\n", t->v_x_mirror[2]);
	printf("v_a_op = %d\n", t->v_x_op[0]);
	printf("v_b_op = %d\n", t->v_x_op[1]);
	printf("v_c_op = %d\n", t->v_x_op[2]);

	bichromatic_decompose(b, r, t, sp, f1, f2);

	return 0;
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2){
	//definitions
	int v_x[3];
	int v_x_next[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	int k;
	for (k = 0; k < 3 && t->vertex_tripod_assign[v_x[k]] != sp; k++){
		//this loop checks who is empty, k == 3 if none
		//if k \in {0, 1, 2} then we're in the first if statement
	}
	//use: k, (k+1)%3, (k+2)%3

	if (k < 3){ //if one of our three legs is non-empty (there can be at most one non-empty leg)
		//here we either have case 5.1 or 5.2
		//there can be at most one non-empty leg, since more empty legs lead to trichromatic subproblems
		//the vertices of the crotch of sp are one of exactly three colours
		if (v_x_next[k] != v_x[(k+1)%3] && v_x_next[k] != v_x[(k+2)%3]){ //if the path up the bfs tree from v_x[k] does NOT lead to v_x[(k+1)%3 or v_x[(k+2)%3]
			//here we have case 5.1
			printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
			if (t->vertex_tripod_assign[v_x_next[k]] == t->vertex_tripod_assign[v_x[(k+1)%3]]){ //if v_x_next[k] is the same colour as v_x[(k+1)%3]
				printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_l[k], v_x_mirror[k]);
				bichromatic_tripod( b, r, t, v_x_l[k], v_x_mirror[k]);
				printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_x_r[k], v_x_mirror[(k+2)%3]);
				trichromatic_tripod(b, r, t, f2, v_x_r[k], v_x_mirror[(k+2)%3]);
			}
			else { //if v_x_next[k] is the same colour as v_x[(k+2)%3]
				printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_r[k], v_x_mirror[(k+2)%3]);
				bichromatic_tripod( b, r, t, v_x_r[k], v_x_mirror[(k+2)%3]);
				printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_x_mirror[k], v_x_l[k]);
				trichromatic_tripod(b, r, t, f2, v_x_mirror[k], v_x_l[k]);
			}
		}
		else {
			//here we have case 5.2
			printf("we have one subproblem for sp %d: trichromatic\n", sp);
			printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_x_mirror[k], v_x_mirror[(k+2)%3], f2);
			trichromatic_tripod(b, r, t, v_x_mirror[k], v_x_mirror[(k+2)%3], f2);
		}
	}
	else { //all legs are empty
		//here we either have case 5.3, 5.4, 5.5, or 5.6 (the empty case)
		//the vertices of the crotch of sp are one of exactly two colours

		int m;
		for (m = 0; m < 3 && t->vertex_tripod_assign[v_x[m]] != t->vertex_tripod_assign[v_x[(m+1)%3]]; m++){
			//this loop checks who is of the same colour
		}
		//use: m, (m+1)%3, (m+2)%3, with m and (m+1)%3 of the same colour

		if (m < 3){ //if m and (m+1)%3 are the same colour, with (m+2)%3 a different colour
			//edge of interest is {v_x[m], v_x[(m+1)%3]}
			//mirror triangle of interest is v_x_mirror[m]
			
			//v_x_mirror[m] will have two of its vertices the same colour as v_x[m] and v_x[(m+1)%3]
			//v_x_mirror[m]'s third vertex is either the same colour as v_x[m] and v_x[(m+1)%3], or is uncoloured
			
			if ((b->bt[v_x[(m+1)%3]] == v_x[m] || b->bt[v_x[m]] == v_x[(m+1)%3]) && (f1 != f2)){
				//if v_x[m] is a bfs parent of v_x[(m+1)%3] or v_x[(m+1)%3] is a bfs parent of v_x[m], then we are in case 5.4
				//problem here. we're not always picking the right mirror triangle
				if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+1)%3]);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+1)%3]);
				}
				else if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+2)%3]);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+2)%3]);
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (b->bt[v_x[(m+1)%3]] != v_x[m] && b->bt[v_x[m]] != v_x[(m+1)%3]){ //are we crossing a coloured edge? if not:
					if (t->face_tripod_assign[v_x_mirror[m]] == -1){ //are we about to select a sp on which to recurse? if not:
						//if f1 == f2 then we are in case 5.5
						printf("we have one subproblem for sp %d: monochromatic\n", sp);
						printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[m]);
						monochromatic_tripod( b, r, t, v_x_mirror[m]);
					}
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (t->face_tripod_assign[v_x_mirror[m]] == -1){ //make sure monochromatic subproblem exists
					printf("we have two subproblems for sp %d: monochromatic + bichromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[m]);
					monochromatic_tripod( b, r, t, v_x_mirror[m]);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//different cases here, depending on specific triangle orientation
				if ((b->bt[v_x[(m+1)%3]] == v_x[(m+2)%3] || b->bt[v_x[(m+2)%3]] == v_x[(m+1)%3]) || t->face_tripod_assign[v_x_mirror[(m+1)%3]] > -1){
					if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+2)%3]);
						bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+2)%3]);
					}
				}
				else {
					if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+1)%3]);
						bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+1)%3]);
					}
				}
			}
		}
		else { // if m == 3
			printf("special case involving the boundary of our planar graph\n"); //not necessarily the boundary...
			printf("all vertices are different colours\n");
			printf("are we done?\n"); //looks like it
		}
	}
}

/******************************************************** monochromatic ********************************************************/
int* monochromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1){
	int sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_x[0]] == sp || t->vertex_tripod_assign[t->v_x[1]] == sp || t->vertex_tripod_assign[t->v_x[2]] == sp){ //if one of our legs is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[sp] == sp){
			printf("sperner triangle already found. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
			exit(0);
		}
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	printf("v_a = %d\n", t->v_x[0]);
	printf("v_b = %d\n", t->v_x[1]);
	printf("v_c = %d\n", t->v_x[2]);
	printf("v_a_next = %d\n", t->v_x_next[0]);
	printf("v_b_next = %d\n", t->v_x_next[1]);
	printf("v_c_next = %d\n", t->v_x_next[2]);
	printf("v_a_mirror = %d\n", t->v_x_mirror[0]);
	printf("v_b_mirror = %d\n", t->v_x_mirror[1]);
	printf("v_c_mirror = %d\n", t->v_x_mirror[2]);

	monochromatic_decompose(b, r, t, sp, f1); //passing sp and f1 is redundant

	return 0;
}

void monochromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1){
	//definitions
	int v_x[3];
	int v_x_next[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	int k;
	for (k = 0; k < 3 && t->vertex_tripod_assign[v_x[k]] != sp; k++){
		//this loop checks who is empty, k == 3 if none
		//if k \in {0, 1, 2} then we're in the first if statement
	}
	//use: k, (k+1)%3, (k+2)%3

	if (k < 3){ //if one of our three legs is non-empty (there can be at most one non-empty leg)
		if (v_x_next[k] != v_x[(k+1)%3] && v_x_next[k] != v_x[(k+2)%3]){ //if the path up the bfs tree from v_x[k] does NOT lead to v_x[(k+1)%3 or v_x[(k+2)%3]
			//here we have case 6.1
			printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
			printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_r[k], v_x_mirror[(k+2)%3]);
			bichromatic_tripod( b, r, t, v_x_r[k], v_x_mirror[(k+2)%3]);
			printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_l[k], v_x_mirror[k]);
			bichromatic_tripod( b, r, t, v_x_l[k], v_x_mirror[k]);
		}
		else {
			//here we have case 6.2
			printf("we have one subproblem for sp %d: bichromatic\n", sp);
			printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_x_mirror[(k+2)%3], v_x_mirror[k]);
			bichromatic_tripod(b, r, t, v_x_mirror[(k+2)%3], v_x_mirror[k]);
		}
	}
	else { //all legs are empty
		//here we either have case 6.3, 6.4, or 6.5 (the empty case)
		//the vertices of sp are one colour

		int m;
		for (m = 0; m < 3 && (b->bt[v_x[(m+1)%3]] != v_x[m] && b->bt[v_x[m]] != v_x[(m+1)%3]); m++){
			//this loop checks whom is a parent of whom
		}
		//use: m, (m+1)%3, (m+2)%3, with v_x[m] and v_x[(m+1)%3] having a child-parent relation

		if (m < 3){ //if m and (m+1)%3 have a child-parent relation
			if (b->bt[v_x[(m+1)%3]] == v_x[m]) { //v_x[m] is a bfs parent of v_x[(m+1)%3]
				if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[(m+1)%3]);
					monochromatic_tripod( b, r, t, v_x_mirror[(m+1)%3]);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //if v_x[(m+1)%3] is a bfs parent of v_x[m]
				//here, if the subproblem is empty, a given edge with purple endpoints is either an edge of the purple tripod, or it is adjacent to a sperner triangle
				if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[(m+2)%3]);
					monochromatic_tripod( b, r, t, v_x_mirror[(m+2)%3]);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else { //if m == 3 (if there is no child-parent relation)
			//here we are in case 6.3

			int p;
			for (p = 0; p < 3 && t->face_tripod_assign[v_x_mirror[p]] == -1; p++){
				//this loop checks which v_x_mirror is a previous sp
			}
			//use: p, (p+1)%3, (p+2)%3, with v_x[p] being the previous sp

			if (p < 3){
				printf("we have two subproblems for sp %d: monochromatic + monochromatic\n", sp);
				//recurse on v_x_mirror[(p+1)%3] and v_x_mirror[(p+2)%3]
				if (t->face_tripod_assign[v_x_mirror[(p+1)%3]] == -1){ //make sure v_x_mirror[(p+1)%3] is not a sp
					printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[(p+1)%3]);
					monochromatic_tripod( b, r, t, v_x_mirror[(p+1)%3]);
				}
				if (t->face_tripod_assign[v_x_mirror[(p+2)%3]] == -1){ //make sure v_x_mirror[(p+2)%3] is not a sp
					printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[(p+2)%3]);
					monochromatic_tripod( b, r, t, v_x_mirror[(p+2)%3]);
				}
			}
			else { //if p == 3
				printf("else statement which we should never reach!!!!\n");
			}
		}
	}
}

void store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp){
	int u;
	for (int i = 0; i < 3; i++){
		u = b->sim[sp][i];
		t->v_x[i] = u;
		t->v_x_next[i] = -1;
		while (t->vertex_tripod_assign[u] == -1){
			t->vertex_tripod_assign[u] = sp;
			t->v_x[i] = u;
			u = b->bt[u];
			t->v_x_next[i] = u;
		}
		t->v_x_mirror[i] = b->tri[sp][i];
	}
}

void tprint(struct tripod_decomposition_struct* t){
	printf("\n");
	printf("v_a = %d\n", t->v_x[0]);
	printf("v_b = %d\n", t->v_x[1]);
	printf("v_c = %d\n", t->v_x[2]);
	printf("v_a_next = %d\n", t->v_x_next[0]);
	printf("v_b_next = %d\n", t->v_x_next[1]);
	printf("v_c_next = %d\n", t->v_x_next[2]);
	printf("v_a_mirror = %d\n", t->v_x_mirror[0]);
	printf("v_b_mirror = %d\n", t->v_x_mirror[1]);
	printf("v_c_mirror = %d\n", t->v_x_mirror[2]);
	printf("v_a_op = %d\n", t->v_x_op[0]);
	printf("v_b_op = %d\n", t->v_x_op[1]);
	printf("v_c_op = %d\n", t->v_x_op[2]);
	printf("\n");
}

void tripod_free(struct tripod_decomposition_struct* t){
	free(t->vertex_tripod_assign);
	free(t->face_tripod_assign);
	free(t->tripod_assign_order);
}

int three_tree_test(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp){
	//sanity check for 3-tree property

	int leg_colours[3];
	for (int i = 0; i < 3; i++){ //identify the endpoint colour of each leg
		(t->v_x_next[i] == -1) ? (leg_colours[i] = t->vertex_tripod_assign[t->v_x[i]]) : (leg_colours[i] = t->vertex_tripod_assign[t->v_x_next[i]]);
	}

	int vertex_index = -1;
	for (int i = 0; i < 3; i++){ //for each leg colour of our sp
		for (int j = 0; j < t->tripod_assign_order_index; j++){ //for each previously assigned tripod
			if (leg_colours[i] == t->tripod_assign_order[j]) vertex_index = j; //record index of parent tripod of that colour
		}
		if (vertex_index == -1) return 0; //if the parent tripod doesn't exist, return false
		else if (vertex_index > t->tripod_assign_order_index) return 0; //if the parent index is larger than sp's index, return false
	}

	return 1;
}

int three_tree_test_pt2(struct bfs_struct* b, struct tripod_decomposition_struct* t){
	return 1;
}
//for each vertex v of sp adjacent to another tripod
//look at the tripod v belongs to
//should be same tripod as v, or one of v's tripod parents (tripods on the boundary when we found the tripod)
//all parents have smaller index
//everyone has at most three parents
//something else to check (with a new function)
//for each vertex make sure the parents are okay, and that a "parent" is not in fact a "grand-parent"
//use the three_tree_test function to keep a list of vertices and their parents (like an ajdacency list but for ancestry)
//then, at the end of the program, when the recursion terminates, call the new function
//this new function is similar to three_tree_test, but it should catch an extra case where v has a parent beyond the boundary of its three defining tripods
//so, for each vertex v, look at the parent of tripod v belongs to
//this should be same tripod as v, or one of v's tripod parents (tripods on the boundary when we found the tripod)
//all parents have smaller index
//everyone has at most three parents
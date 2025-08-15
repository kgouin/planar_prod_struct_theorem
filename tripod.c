#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void init(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, char adjacencies[], char simplices[], char triangles[]){
	//initialize bfs_struct & rmq_struct
	BFS_init(b, adjacencies, simplices, triangles);
	BFS(b);
	LCA_init(r, b->ct, (2*(b->f))-1);

	//initialize tripod_decomposition_struct
	t->vertex_tripod_assign = malloc((b->v)*sizeof(int));
	//label the three vertices incident to face 0 as belonging to three exterior tripods
	t->vertex_tripod_assign[b->sim[0][0]] = (b->f);
	t->vertex_tripod_assign[b->sim[0][1]] = (b->f)+1;
	t->vertex_tripod_assign[b->sim[0][2]] = (b->f)+2;
	for (int i = 3; i < (b->v); i++){
		t->vertex_tripod_assign[i] = -1;
	}
	t->face_tripod_assign = malloc((b->f)*sizeof(int));
	t->face_tripod_assign[0] = 0; //this is a special case, since 0 is the outer face
	for (int i = 1; i < (b->f); i++){
		t->face_tripod_assign[i] = -1;
	}
	t->tripod_adjacency_list = malloc((b->f+3)*sizeof(int*));
	for (int i = 0; i < (b->f+3); i++){
		t->tripod_adjacency_list[i] = malloc(3*sizeof(int));
		for (int j = 0; j < 3; j++){
			t->tripod_adjacency_list[i][j] = -1;
		}
	}//we will eventually fill t->tripod_adjacency_list with tripods contributing at least one vertex in H. the other tripods will have values of -1
	t->tripod_adjacency_list[(b->f+1)][0] = (b->f);
	t->tripod_adjacency_list[(b->f+2)][0] = (b->f);
	t->tripod_adjacency_list[(b->f+2)][1] = (b->f+1);
	t->tripod_assign_order = malloc((b->f+3)*sizeof(int));
	//label the first three entries as the three exterior tripods that constitute our base case
	t->tripod_assign_order[0] = (b->f);
	t->tripod_assign_order[1] = (b->f)+1;
	t->tripod_assign_order[2] = (b->f)+2;
	for (int i = 3; i < (b->f+3); i++){
		t->tripod_assign_order[i] = -1;
	}
	t->tripod_assign_order_index = 3;
	//label portals in clockwise, not counterclockwise order, since we are dealing with the outer face
	//label portals according to their indices in simplices
	t->portals[0][0] = b->tri[0][0];
	int j;
	for (j = 0; j < 3 && b->sim[t->portals[0][0]][j] != 1; j++){
		//find the index of vertex 1 in simplices at index [outer face]
	}
	//from there, store indices in t->portals array
	t->portals[0][1] = j;
	t->portals[1][0] = b->tri[0][1];
	t->portals[1][1] = (j+1)%3;
	t->portals[2][0] = b->tri[0][2];
	t->portals[2][1] = (j+2)%3;
}

void decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){
	//start decomposition with the three triangles adjacent to outer face
	trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1]); //these are in counterclockwise order
}

/******************************************************** less intrusive portal helpers ********************************************************/

void less_intrusive_portal_helper_r_l(struct bfs_struct* b, struct tripod_decomposition_struct* t, int r_l_triangle, int v, int portal_index){
	//as input to this function, for the left triangle we want to pass v_x and for the right triangle we want to pass v_x_next
	t->portals[portal_index][0] = r_l_triangle;
	int j;
	for (j = 0; j < 3 && b->sim[t->portals[portal_index][0]][j] != v; j++){
		//find the index of vertex v in simplices at index [t->portals[portal_index][0]]
	}
	t->portals[portal_index][1] = j;
}

void less_intrusive_portal_helper_op(struct bfs_struct* b, struct tripod_decomposition_struct* t, int portals[3][2], int op_triangle, int portal_index){
	t->portals[portal_index][0] = op_triangle;
	for (int j = 0; j < 3; j++){ //look through previous portals to find portal pair [op_triangle, index of v]. we want this to be our current portal
		if (portals[j][0] == t->portals[portal_index][0]){
			t->portals[portal_index][1] = portals[j][1];
			break;
		}
	}
}

void less_intrusive_portal_helper_mirror(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int mirror_triangle, int portal_index){
	int p;
	t->portals[portal_index][0] = mirror_triangle;
	for (int j = 0; j < 3; j++){ //look through b->tri at index sp to find index j of mirror_triangle in that array. v will be b->sim[sp][j]
		if (b->tri[sp][j] == t->portals[portal_index][0]){
			p = b->sim[sp][(j+1)%3]; //this is possible because of the way our data is arranged
			break;
		}
	}

	int j;
	for (j = 0; j < 3 && b->sim[t->portals[portal_index][0]][j] != p; j++){
		//find the index of vertex v in simplices at index [t->portals[portal_index][0]]
	}
	t->portals[portal_index][1] = j;
}

void less_intrusive_portal_helper_empty(struct tripod_decomposition_struct* t, int portal_index){
	t->portals[portal_index][0] = -1;
	t->portals[portal_index][1] = -1;
}

/******************************************************** trichromatic ********************************************************/

void trichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3){
	//sperner triangle identification
	int sp;
	if (f1 == f2 && f1 == f3){
		sp = f1;
		t->face_tripod_assign[sp] = sp; //keep track of sperner triangles
		//don't add tripod to t->tripod_assign_order
		return;
	}
	else if (f1 != f2 && f1 != f3 && f2 != f3){
		//find sperner triangle
		if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3) && LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		else {
			if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3)) sp = LCA_query(r, f2, f3);
			else if (LCA_query(r, f1, f2) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f3);
			else if (LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
			else exit(0);
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

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	trichromatic_orient_subproblems(b, t, sp, f1, f2, f3);

	trichromatic_decompose(b, r, t, sp, f1, f2, f3);

	return;
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

	int double_match = 0;
	for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
		if (t->vertex_tripod_assign[b->sim[f[k]][i]] == leg_colours[0]){ //this can be simpler with portal edges
			for (int j = 0; j < 3; j++){
				if (t->vertex_tripod_assign[b->sim[f[k]][j]] == leg_colours[1]){
					double_match = 1;
				}
			}
		}
	}
	if (double_match){
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
			t->v_x_op[0] = f[(k+2)%3];
			t->v_x_op[1] = f[k];
			t->v_x_op[2] = f[(k+1)%3];
		}
		else { //no double match. f must match with leg_colour_c and leg_colour_a
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
	int portals[3][2];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_op[i] = t->v_x_op[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
		portals[i][0] = t->portals[i][0];
		portals[i][1] = t->portals[i][1];
	}

	for (int i = 0; i < 3; i++){
		if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is non-empty && leg (i+1)%3 is non-empty
			less_intrusive_portal_helper_op(b, t, portals, v_x_op[i], 0);
			less_intrusive_portal_helper_r_l(b, t, v_x_r[(i+1)%3], v_x_next[(i+1)%3], 1);
			less_intrusive_portal_helper_r_l(b, t, v_x_l[i], v_x[i], 2);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is non-empty && leg (i+1)%3 is empty
			less_intrusive_portal_helper_op(b, t, portals, v_x_op[i], 0);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[i], 1);
			less_intrusive_portal_helper_r_l(b, t, v_x_l[i], v_x[i], 2);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_mirror[i], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is empty && leg (i+1)%3 is non-empty
			less_intrusive_portal_helper_op(b, t, portals, v_x_op[i], 0);
			less_intrusive_portal_helper_r_l(b, t, v_x_r[(i+1)%3], v_x_next[(i+1)%3], 1);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[i], 2);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_mirror[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is empty && leg (i+1)%3 is empty
			if (t->face_tripod_assign[v_x_op[i]] == -1){
				less_intrusive_portal_helper_op(b, t, portals, v_x_op[i], 0);
				less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[i], 1);
				less_intrusive_portal_helper_empty(t, 2);
				bichromatic_tripod( b, r, t, v_x_op[i], v_x_mirror[i]);
			}
		}
	}
}

/******************************************************** bichromatic ********************************************************/

void bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2){
	int sp;

	sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_x[0]] == sp || t->vertex_tripod_assign[t->v_x[1]] == sp || t->vertex_tripod_assign[t->v_x[2]] == sp){ //if one of our legs is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	bichromatic_decompose(b, r, t, sp, f1, f2);
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2){
	//definitions
	int v_x[3];
	int v_x_next[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	int portals[3][2];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
		portals[i][0] = t->portals[i][0];
		portals[i][1] = t->portals[i][1];
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
			if (t->vertex_tripod_assign[v_x_next[k]] == t->vertex_tripod_assign[v_x[(k+1)%3]]){ //if v_x_next[k] is the same colour as v_x[(k+1)%3]
				less_intrusive_portal_helper_r_l(b, t, v_x_l[k], v_x[k], 0);
				less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[k], 1);
				less_intrusive_portal_helper_empty(t, 2);
				bichromatic_tripod( b, r, t, v_x_l[k], v_x_mirror[k]);

				less_intrusive_portal_helper_op(b, t, portals, f2, 0);
				less_intrusive_portal_helper_r_l(b, t, v_x_r[k], v_x_next[k], 1);
				less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(k+2)%3], 2);
				trichromatic_tripod(b, r, t, f2, v_x_r[k], v_x_mirror[(k+2)%3]);
			}
			else { //if v_x_next[k] is the same colour as v_x[(k+2)%3]
				less_intrusive_portal_helper_r_l(b, t, v_x_r[k], v_x_next[k], 0);
				less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(k+2)%3], 1);
				less_intrusive_portal_helper_empty(t, 2);
				bichromatic_tripod( b, r, t, v_x_r[k], v_x_mirror[(k+2)%3]);

				less_intrusive_portal_helper_op(b, t, portals, f2, 0);
				less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[k], 1);
				less_intrusive_portal_helper_r_l(b, t, v_x_l[k], v_x[k], 2);
				trichromatic_tripod(b, r, t, f2, v_x_mirror[k], v_x_l[k]);
			}
		}
		else {
			//here we have case 5.2
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[k], 0);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(k+2)%3], 1);
			less_intrusive_portal_helper_op(b, t, portals, f2, 2);
			trichromatic_tripod(b, r, t, v_x_mirror[k], v_x_mirror[(k+2)%3], f2);
		}
	}
	else { //all legs are empty
		//here we either have case 5.3, 5.4, 5.5, or 5.6 (the empty case)
		//the vertices of the crotch of sp are one of exactly two colours

		int m = -1;
		for (int j = 0; j < 3; j++){
			if (portals[j][0] == sp){
				m = portals[j][1];
				break;
			}
		}

		if (m < 3){ //this should always be the case

			//v_x_mirror[m] will have two of its vertices the same colour as v_x[m] and v_x[(m+1)%3]
			//v_x_mirror[m]'s third vertex is either the same colour as v_x[m] and v_x[(m+1)%3], or is uncoloured

			if (((b->bt[v_x[(m+1)%3]] == v_x[(m+2)%3] || b->bt[v_x[(m+2)%3]] == v_x[(m+1)%3]) || (b->bt[v_x[(m+2)%3]] == v_x[m] || b->bt[v_x[m]] == v_x[(m+2)%3])) && (f1 != f2)){
				//if there is a parent-child relation, then we are in the first case of 5.4, where we have ancestry
				//although this second case of 5.4 is really just a case 5.3, with no monochromatic problem on which to recurse
				//second case of 5.4 is dealt with later

				//use portal edges here

				int ell;

				if (b->bt[v_x[(m+1)%3]] == v_x[(m+2)%3] || b->bt[v_x[(m+2)%3]] == v_x[(m+1)%3]){
					ell = (m+2)%3;
				}
				else if (b->bt[v_x[(m+2)%3]] == v_x[m] || b->bt[v_x[m]] == v_x[(m+2)%3]){
					ell = (m+1)%3;
				}
				else exit(0);

				if (t->face_tripod_assign[v_x_mirror[ell]] == -1){ //if face to be selected doesn't already belong to a tripod
					less_intrusive_portal_helper_op(b, t, portals, f2, 0);
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[ell], 1);
					less_intrusive_portal_helper_empty(t, 2);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[ell]);
				}
				//if face to be selected does already belong to a tripod, then there is nothing left on which to recurse
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				//if f1 == f2 then we are in case 5.5
				//here we have two valid values of m for the same triangle
				//find the second value of m:
				int second_m = -1;
				for (int j = 0; j < 3; j++){
					if (portals[j][0] == sp && portals[j][1] != m){
						second_m = portals[j][1];
						break;
					}
				}

				int ell = 3-(m+second_m);

				if (t->face_tripod_assign[v_x_mirror[ell]] == -1 && (b->bt[v_x[ell]] != v_x[(ell+1)%3] && b->bt[v_x[(ell+1)%3]] != v_x[ell])){
					//if face to be selected doesn't already belong to a tripod AND there is no parent-child relation between vertices of edge whose other side we're about to recurse on
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[ell], 0);
					less_intrusive_portal_helper_empty(t, 1);
					less_intrusive_portal_helper_empty(t, 2);
					monochromatic_tripod( b, r, t, v_x_mirror[ell]);
				}
				//if face to be selected does already belong to a tripod, then there is nothing left on which to recurse
			}
			else {
				//we are in case 5.3
				//or we are in the second case of 5.4 where our new sp is adjacent to a previous sp, and by extension has no parent-child relation
				//as noted above, this second case of 5.4 is really just a case 5.3, with no monochromatic problem on which to recurse
				//for any coloured edge, either it's a tree edge or there's a sp on the other side

				if (t->vertex_tripod_assign[b->sim[sp][(m+1)%3]] == t->vertex_tripod_assign[b->sim[sp][(m+2)%3]]){ //this is our monochromatic edge
					if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){ //if a monochromatic subproblem exists
						less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+1)%3], 0);
						less_intrusive_portal_helper_empty(t, 1);
						less_intrusive_portal_helper_empty(t, 2);
						monochromatic_tripod( b, r, t, v_x_mirror[(m+1)%3]);
					}
					less_intrusive_portal_helper_op(b, t, portals, f2, 0);
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+2)%3], 1);
					less_intrusive_portal_helper_empty(t, 2);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+2)%3]);
				}
				else if (t->vertex_tripod_assign[b->sim[sp][(m+2)%3]] == t->vertex_tripod_assign[b->sim[sp][m]]){ //this is our monochromatic edge
					if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){ //if a monochromatic subproblem exists
						less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+2)%3], 0);
						less_intrusive_portal_helper_empty(t, 1);
						less_intrusive_portal_helper_empty(t, 2);
						monochromatic_tripod( b, r, t, v_x_mirror[(m+2)%3]);
					}
					less_intrusive_portal_helper_op(b, t, portals, f2, 0);
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+1)%3], 1);
					less_intrusive_portal_helper_empty(t, 2);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+1)%3]);
				}
				else exit(0);
			}
		}
		else exit(0); //if m == 3
	}
}

/******************************************************** monochromatic ********************************************************/
void monochromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1){
	int sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_x[0]] == sp || t->vertex_tripod_assign[t->v_x[1]] == sp || t->vertex_tripod_assign[t->v_x[2]] == sp){ //if one of our legs is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	monochromatic_decompose(b, r, t, sp);
}

void monochromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp){
	//definitions
	int v_x[3];
	int v_x_next[3];
	int v_x_mirror[3];
	int v_x_l[3];
	int v_x_r[3];
	int portals[3][2];
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
		portals[i][0] = t->portals[i][0];
		portals[i][1] = t->portals[i][1];
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
			less_intrusive_portal_helper_r_l(b, t, v_x_r[k], v_x_next[k], 0);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(k+2)%3], 1);
			less_intrusive_portal_helper_empty(t, 2);
			bichromatic_tripod( b, r, t, v_x_r[k], v_x_mirror[(k+2)%3]);

			less_intrusive_portal_helper_r_l(b, t, v_x_l[k], v_x[k], 0);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[k], 1);
			less_intrusive_portal_helper_empty(t, 2);
			bichromatic_tripod( b, r, t, v_x_l[k], v_x_mirror[k]);
		}
		else {
			//here we have case 6.2
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(k+2)%3], 0);
			less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[k], 1);
			less_intrusive_portal_helper_empty(t, 2);
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
			if (b->bt[v_x[(m+1)%3]] == v_x[m]){ //v_x[m] is a bfs parent of v_x[(m+1)%3]
				if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+1)%3], 0);
					less_intrusive_portal_helper_empty(t, 1);
					less_intrusive_portal_helper_empty(t, 2);
					monochromatic_tripod( b, r, t, v_x_mirror[(m+1)%3]);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //if v_x[(m+1)%3] is a bfs parent of v_x[m]
				//here, if the subproblem is empty, a given edge with purple endpoints is either an edge of the purple tripod, or it is adjacent to a sperner triangle
				if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(m+2)%3], 0);
					less_intrusive_portal_helper_empty(t, 1);
					less_intrusive_portal_helper_empty(t, 2);
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
				//recurse on v_x_mirror[(p+1)%3] and v_x_mirror[(p+2)%3]
				if (t->face_tripod_assign[v_x_mirror[(p+1)%3]] == -1){ //make sure v_x_mirror[(p+1)%3] is not a sp
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(p+1)%3], 0);
					less_intrusive_portal_helper_empty(t, 1);
					less_intrusive_portal_helper_empty(t, 2);
					monochromatic_tripod( b, r, t, v_x_mirror[(p+1)%3]);
				}
				if (t->face_tripod_assign[v_x_mirror[(p+2)%3]] == -1){ //make sure v_x_mirror[(p+2)%3] is not a sp
					less_intrusive_portal_helper_mirror(b, t, sp, v_x_mirror[(p+2)%3], 0);
					less_intrusive_portal_helper_empty(t, 1);
					less_intrusive_portal_helper_empty(t, 2);
					monochromatic_tripod( b, r, t, v_x_mirror[(p+2)%3]);
				}
			}
			else exit(0); //if p == 3
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
		t->tripod_adjacency_list[sp][i] = t->vertex_tripod_assign[u]; //t->vertex_tripod_assign[u] is colour of parent 1 of possibly 3
		t->v_x_mirror[i] = b->tri[sp][i];
	}
}

void tripod_free(struct bfs_struct* b, struct tripod_decomposition_struct* t){
	free(t->vertex_tripod_assign);
	free(t->face_tripod_assign);
	free(t->tripod_assign_order);
	for (int i = 0; i < (b->f+3); i++){
		free(t->tripod_adjacency_list[i]);
	}
	free(t->tripod_adjacency_list);
}

int three_tree_test(struct bfs_struct* b, struct tripod_decomposition_struct* t){
	//second sanity check for 3-tree property
	//for each edge, do its endpoints belong to the same tripod, or is one the parent of the other: then all good

	int pass = 0;

	for (int i = 0; i < (b->v); i++){ //go through all vertices
		for (int j = 0; j < (b->n)[i]; j++){ //look at adjacency list of each vertex, using (b->n)[i] as the number of neighbours of vertex i
			//use b->al[i][j] as i's jth neighbour
			//use t->vertex_tripod_assign[i] as i's tripod colour
			//use t->vertex_tripod_assign[b->al[i][j]] as i's jth neighbour's tripod colour
			if (t->vertex_tripod_assign[i] == t->vertex_tripod_assign[b->al[i][j]]) pass = 1; //if i's tripod is the same as j's tripod, then ok
			for (int k = 0; k < 3; k++){
				if (t->tripod_adjacency_list[t->vertex_tripod_assign[b->al[i][j]]][k] == t->vertex_tripod_assign[i]) pass = 1; //if i's tripod is in j's tripod adjacency list, then ok
				if (t->tripod_adjacency_list[t->vertex_tripod_assign[i]][k] == t->vertex_tripod_assign[b->al[i][j]]) pass = 1; //if j's tripod is in i's tripod adjacency list, then ok
			}
			if (!pass) return 0;
			pass = 0;
		}
	}
	return 1;
}

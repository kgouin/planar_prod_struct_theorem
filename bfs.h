struct bfs_struct{
	int v; //the number of vertices in a triangulation
	int** a; //the adjacency list of a triangulation
	int* n; //the number of neighbours for each vertex in a triangulation
	int* r; //the roots of the triangulation
	int* bfs; //the breadth-first-search array
};

void BFS_init(struct bfs_struct*);
int* BFS(struct bfs_struct*);
void BFS_free(struct bfs_struct*);

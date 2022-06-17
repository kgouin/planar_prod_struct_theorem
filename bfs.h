struct bfs_struct{
	int v; //the number of vertices in a triangulation
	int** a; //the adjacency list of a triangulation
	int* n; //the number of neighbours for each vertex in a triangulation
	int* p; //the parent array, which will act as the seen array as well
};

void BFS_init(struct bfs_struct*);
int* BFS(struct bfs_struct*);
void BFS_free(struct bfs_struct*);

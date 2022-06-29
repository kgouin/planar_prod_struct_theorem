struct bfs_struct{
	int v; //the number of vertices in a triangulation
	int f; //the number of faces in a triangulation
	int* n; //the number of neighbours for each vertex in a triangulation
	int** al; //the adjacency list of a triangulation (adjacencies.txt information)
	int** il; //the incidence list of a triangulation (adjacencies.txt information)
	int** sim; //simplices.txt information
	int** tri; //triangles.txt information
	int* bt; //the breadth-first search spanning tree representation containing parents of visited vertices
	int* pin; //the breadth-first search spanning tree representation containing the indices of the parents of visited vertices in al
	int** ct; //the cotree representation
};

void BFS_init(struct bfs_struct*);
int* BFS(struct bfs_struct*);
void BFS_free(struct bfs_struct*);

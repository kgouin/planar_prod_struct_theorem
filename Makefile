
CFLAGS=-g -Wall   # for debugging
CFLAGS=-g -Wall -pg  # for debugging and profiling
CFLAGS=-O -Wall  # fastest code
all : rmq_exec lca_exec bfs_exec tripod_exec clean

rmq_exec : rmq.h rmq rmq_test.c
	gcc $(CFLAGS) -o rmq_exec rmq.c rmq_test.c -lm

rmq : rmq.h
	gcc $(CFLAGS) -c rmq.c

lca_exec : lca.h rmq.h lca lca_test.c
	gcc $(CFLAGS) -o lca_exec rmq.c lca.c lca_test.c -lm

lca : rmq.h lca.h
	gcc $(CFLAGS) -c lca.c

bfs_exec : bfs.h bfs bfs_test.c
	gcc $(CFLAGS) -o bfs_exec bfs.c bfs_test.c -lm

bfs : bfs.h
	gcc -c bfs.c

tripod : rmq.h lca.h bfs.h tripod.h
	gcc $(CFLAGS) -c tripod.c

tripod_exec : rmq.h lca.h bfs.h tripod.h tripod tripod_test.c
	gcc $(CFLAGS) -o tripod_exec rmq.c lca.c bfs.c tripod.c tripod_test.c -lm

profile.txt : rmq.h rmq rmq_test.c
	gcc $(CFLAGS) -o rmq_prof rmq.c rmq_test.c -lm
	./rmq_prof
	gprof -l rmq_prof gmon.out > profile.txt
	echo "See profile.txt for line-level profiling information"

clean :
	rm -f *.o

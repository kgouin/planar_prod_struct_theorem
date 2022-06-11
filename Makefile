all : rmq_exec lca_exec clean

rmq_exec : rmq.h rmq rmq_test.c
	gcc -g -o rmq_exec rmq.c rmq_test.c -lm

rmq : rmq.h
	gcc -c rmq.c

lca_exec : lca.h rmq.h lca lca_test.c
	gcc -g -o lca_exec rmq.c lca.c lca_test.c -lm

lca : rmq.h lca.h
	gcc -c lca.c

profile.txt : rmq.h rmq rmq_test.c
	gcc -g -pg -o rmq_prof rmq.c rmq_test.c -lm
	./rmq_prof
	gprof -l rmq_prof gmon.out > profile.txt
	echo "See profile.txt for line-level profiling information"

clean :
	rm -f *.o

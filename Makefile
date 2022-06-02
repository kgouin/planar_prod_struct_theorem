


lca : lca.c
	gcc -g -o lca lca.c -lm



profile.txt : lca.c
	gcc -g -pg -o lca_prof lca.c -lm
	./lca_prof
	gprof -l lca_prof gmon.out > profile.txt
	echo "See profile.txt for line-level profiling information"

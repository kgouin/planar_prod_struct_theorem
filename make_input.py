#!/usr/bin/python3
import sys
import time
import random
import itertools
import collections
import matplotlib.pyplot as plt
import scipy.spatial


def triangulate(points):
    n = len(points)
    dt = scipy.spatial.Delaunay(points)

    assert(dt.npoints == n)
    assert(len(dt.convex_hull) == 3)
    assert(dt.nsimplex == 2*n - 5)
    left_of = [collections.defaultdict(int) for _ in range(n)]
    succ = [dict() for _ in range(n)]
    simplices = [[-1]*3] + [list(t) for t in dt.simplices]
    for k in range(1, len(simplices)):
        t = simplices[k]
        for i in range(3):
            succ[t[i]][t[(i+1)%3]] = t[(i+2)%3]
            left_of[t[i]][t[(i+1)%3]] = k

    triangles = [[-1]*3 for _ in range(len(simplices))]
    for k in range(1, len(simplices)):
        t = simplices[k]
        for i in range(3):
            if t[i] in left_of[t[(i+1)%3]]:
                triangles[k][i] = left_of[t[(i+1)%3]][t[i]]

    # don't forget the outer face
    of = list(set(itertools.chain.from_iterable(dt.convex_hull)))

    if of[1] in left_of[of[0]]:
        of = of[::-1]
    simplices[0] = of[:]
    for i in range(3):
        t2 = left_of[of[(i+1)%3]][of[i]]
        #simplices[0][i] = succ[of[(i+1)%3]][of[i]]
        succ[of[i]][of[(i+1)%3]] = of[(i+2)%3]
        triangles[0][i] = t2
        triangles[t2][triangles[t2].index(-1)] = 0

    # Each triangle stores indices of three adjacent triangles
    print("Writing triangles.txt")
    with open("triangles.txt", "w") as fp:
        fp.write("{}\n".format(len(triangles)))
        for t in triangles:
            fp.write("{} {} {}\n".format(t[0], t[1], t[2]))

    # Each triangle stores indices of three incident vertices
    print("Writing simplices.txt")
    with open("simplices.txt", "w") as fp:
        fp.write("{}\n".format(len(simplices)))
        for t in simplices:
            fp.write("{} {} {}\n".format(t[0], t[1], t[2]))

    # Each vertex stores its adjacency list
    print("Writing adjacencies.txt")
    al = succ2al(succ)
    with open("adjacencies.txt", "w") as fp:
        fp.write("{}\n".format(n))
        print(al)
        for i in range(n):
            # length of adjacency list
            fp.write("{} ".format(len(al[i])))

            # adjacency list (vertices)
            fp.write("{} ".format(" ".join([str(j) for j in al[i]])))

            # incidence list (triangles)
            print(al[i])
            fp.write("{}\n".format(" ".join([str(left_of[i][j]) for j in al[i]])))



    # Each vertex stores the coordinates of the point that defines it
    print("Writing points.txt")
    with open("points.txt", "w") as fp:
        fp.write("{}\n".format(n))
        for p in points:
            fp.write("{} {}\n".format(p[0], p[1]))






######################################################################
# Boring routines to build a "random" triangulation
######################################################################
def make_triangulation(n, data_type):
    print("Generating points")
    if data_type == 0:
        # Use a set of n-3 random points
        points = [(-1.5,-1.5), (-1.5,3), (3,-1.5)] \
                 + [random_point() for _ in range(n-3)]
    elif data_type == 1:
        # Use a set of n-3 collinear points
        points = [(-1.5,-1.5), (-1.5,3), (3,-1.5)] \
                 + [(-1 + i/(n-3), -1 + i/(n-3)) for i in range(n-3)]
    elif data_type == 2:
        points = [(0, 0), (1,1), (1,0)] \
                 + [(random.random(), random.random()) for _ in range(n-3)]
        for i in range(n):
            (x, y) = points[i]
            if x < y:
                points[i] = (y, x)
    else:
        raise ValueError("Invalid argument for data_type")

    n = len(points)
    # random.shuffle(points)

    print("Computing Delaunay triangulation")
    triangulate(points)


""" Generate a random point in the unit circle """
def random_point():
    while 1 < 2:
        x = 2*random.random()-1
        y = 2*random.random()-1
        if x**2 + y**2 < 1:
            return (x, y)

""" Convert a triangle-based adjacency representation into an adjacency-list representation """
def succ2al(succ):
    al = list()
    for sd in succ:
        al.append(list())
        v0 = next(iter(sd))
        v = v0
        while True: # emulating do ... while v != v0
            al[-1].append(v)
            v = sd[v]
            if v == v0: break
    return al

def usage():
    print("Computes a tripod decomposition of a Delaunay triangulation")
    print("Usage: {} [-h] [-c] [-r] [-y] [-w] [-b] [-nv] <n>".format(sys.argv[0]))
    print("  -h show this message")
    print("  -c use collinear points")
    print("  -y use random points in triangle")
    print("  -r use random points in disk (default)")
    print("  -w use O(n log n) time algorithm (default)")
    print("  -b use O(n^2) time algorithm (usually faster)")
    print("  -nv don't verify correctness of results")
    print("  <n> the number of points to use (default = 10)")

if __name__ == "__main__":
    n = 0
    data_type = 0
    worst_case = True
    verify = True
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-r':
            data_type = 0   # random
        elif arg == '-c':
            data_type = 1   # collinear
        elif arg == '-y':
            data_type = 2   # random in triangle (like rbox y)
        elif arg == '-w':
            worst_case = True
        elif arg == '-b':
            worst_case = False
        elif arg == '-nv':
            verify = False
        else:
            n = int(arg)

    if n <= 0:
        usage()
        sys.exit(-1)

    s = ["random", "collinear", "uniform"][data_type]
    print("Generating {} point set of size {}".format(s, n))
    make_triangulation(n, data_type)

    # n = len(succ)
    # m = sum([len(x) for x in succ]) // 2
    # print("n = ", n, " m = ", m)
    # assert(m == 3*n - 6)
    #
    # if n > 500:
    #     print("Not displaying results since n = {} > 500".format(n))
    #     sys.exit(0)
    #
    # # Draw graph
    # for v in range(len(succ)):
    #     for w in succ[v]:
    #         plt.plot([points[v][0], points[w][0]], [points[v][1], points[w][1]], color='gray', lw=0.2)
    #
    # for v in range(n):
    #     plt.plot(points[v][0], points[v][1], color="red", lw=1, marker='o',
    #              markersize=min(8,180/n))
    #
    #
    # cmap = ['red', 'darkgreen', 'blue', 'orange', 'ghostwhite']
    # fmap = ['mistyrose', 'lightgreen', 'lightblue', 'moccasin', 'ghostwhite']
    #
    # # Draw tripods
    # tripod_colours = tp.colour_tripods()
    # for i in range(1, len(tp.tripods)):
    #     tripod = tp.tripods[i]
    #     c = tripod_colours[i]
    #     # Draw legs
    #     for path in tripod:
    #         x = [points[v][0] for v in path]
    #         y = [points[v][1] for v in path]
    #         plt.plot(x, y, color=cmap[c], lw=2)
    #     tau = [tripod[i][0] for i in range(3)]
    #     # Draw and label Sperner triangle
    #     x = [points[v][0] for v in tau]
    #     y = [points[v][1] for v in tau]
    #     if n <= 100:
    #         plt.fill(x, y, facecolor=fmap[c], lw=0)
    #     x = sum(x)/3
    #     y = sum(y)/3
    #     if n < 250:
    #         plt.text(x, y, str(i), horizontalalignment='center',
    #                  verticalalignment='center', fontsize=min(10,500/n))
    #
    #     tau2 = sum([tripod[j][:-1][:1] for j in range(3)], [])
    #     if tau2:
    #         tau2.append(tau2[0])
    #         x = [points[v][0] for v in tau2]
    #         y = [points[v][1] for v in tau2]
    #         plt.plot(x, y, color=cmap[c], lw=2)
    #
    # for v in range(n):
    #     t = tp.tripod_map[v][0]
    #     c = cmap[tripod_colours[t]]
    #     plt.plot(points[v][0], points[v][1], color=c, lw=1, marker='o',
    #              markersize=min(8,400/n))
    #
    # plt.axis('off')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()

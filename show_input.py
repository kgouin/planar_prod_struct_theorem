#!/usr/bin/python3
import sys
import time
import random
import itertools
import matplotlib.pyplot as plt
import scipy.spatial


if __name__ == "__main__":
    triangles = list()
    with open("simplices.txt", "r") as fp:
        fp.readline()
        for line in fp:
            triangles.append([int(x) for x in line.split()])
    #print(triangles)
    with open("points.txt", "r") as fp:
        fp.readline()
        points = list()
        for line in fp:
            points.append( tuple([float(x) for x in line.split()]) )
    #print(points)
    bfs = list()
    with open("bfs.txt", "r") as fp:
        fp.readline()
        for line in fp:
            bfs.append(int(line))
    #print(bfs)

    n = len(points)
    #print(n, len(triangles))

    # for v in range(len(succ)):
    #     for w in succ[v]:
    #         plt.plot([points[v][0], points[w][0]], [points[v][1], points[w][1]], color='gray', lw=0.2)

    for t in range(len(triangles)):
        tau = triangles[t]
        #print(tau)
        # Draw and label Sperner triangle
        #for v in tau:
            #print(points[v])
        x = [points[v][0] for v in tau]
        y = [points[v][1] for v in tau]
        # if n <= 100:
        #     plt.fill(x, y, facecolor=fmap[c], lw=0)
        plt.plot(x + [x[0]], y + [y[0]], color='0.8', lw=1)
        x = sum(x)/3
        y = sum(y)/3
        if t > 0 and n < 250:
            plt.text(x, y, str(t), horizontalalignment='center', verticalalignment='center', fontsize=min(8,1000/n), color='tab:orange')

    for v in range(len(points)):
        x,y = points[v]
        first_coord = [x, y]
        # plt.plot(x,y, color="red", lw=1, marker='o',
        #          markersize=min(8,180/n))
        plt.text(x, y, str(v), horizontalalignment='center', verticalalignment='center', fontsize=min(8,500/n))
        for i in range(len(points)):
            if (bfs[i] == v):
                x2, y2 = points[i]
                second_coord = [x2, y2]
                #print("draw from coord " + str(first_coord) + " to coord " + str(second_coord))
                plt.plot([first_coord[0], x2], [first_coord[1], y2], color='tab:purple', lw=1.5)

    plt.axis('off')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

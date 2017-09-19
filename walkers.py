import sys, os
import numpy as np
import random
from math import exp
import matplotlib.pyplot as plt
import time
import cProfile
import re



random.seed(1002)

n=200
n_bugs = 20
ising_array = [[-1 if (random.random() > 0.5) else 1  for i in range(n)] for j in range(n)]
bug_coord = [[random.randint(0,n-1) for i in range(n_bugs)] for j in range(2)]
J = 1
H = 0.0
beta = 1

def randomize():
    global ising_array
    for i in range(n):
        for j in range(n):
            if random.random() < 0.5:
                ising_array[i,j] *= -1


def calc_total_energy(J):
    global ising_array
    energy = 0
    for i in range(n):
        for j in range(n):
            energy += -J*ising_array[i][(j+1)%n] - J*ising_array[(i+1)%n][j]# + H*ising_array[i][j]
    return energy


def calc_del_e(spin, J):
    global ising_array
    i = spin[0]
    j = spin[1]
    neigh = [ising_array[(i+1)%n][j], ising_array[i][(j+1)%n], ising_array[(i-1)%n][j], ising_array[i][(j-1)%n]]
    #when ising array is a numpy array neigh = ising_array[[(i+1)%n,i,(i-1)%n,i],[j,(j+1)%n,j,(j-1)%n]]
    sigma_spin = ising_array[i][j]
    old_e = -J*sum([sigma_spin*i for i in neigh])# + H*sigma_spin
    new_e = -J*sum([-sigma_spin*i for i in neigh])# - H*sigma_spin
    return (new_e - old_e)





tot_en = calc_total_energy(J)
sweeps = 500
n_bug_steps = 50
alive = [1 for i in range(n_bugs)]
for s in range(sweeps):
    #update ising base
    for sw in range(n*n):
        spin = [random.randint(0,n-1),random.randint(0,n-1)]
        del_e = calc_del_e(spin, J)
        if(exp(-beta*del_e) > random.random()):
            ising_array[spin[0]][spin[1]] *= -1
            tot_en += del_e
            #move bugs
    food_consumed = [0 for i in range(n_bugs)]
    for bsw in range(n_bug_steps):
        for cbug in range(n_bugs):
            if alive[cbug] == 1:
                #move random step
                if random.random() > 0.5:
                    if random.random() > 0.5:
                        bug_coord[0][cbug] = (bug_coord[0][cbug] + 1)%n
                    else:
                        bug_coord[0][cbug] = (bug_coord[0][cbug] - 1)%n
                else:
                    if random.random() > 0.5:
                        bug_coord[1][cbug] = (bug_coord[1][cbug] + 1)%n
                    else:
                        bug_coord[1][cbug] = (bug_coord[1][cbug] - 1)%n
                #eat up!
                if ising_array[bug_coord[0][cbug]][bug_coord[1][cbug]] == 1:
                    food_consumed[cbug] += 1
                    ising_array[bug_coord[0][cbug]][bug_coord[1][cbug]] = -1
        

    print(food_consumed)
    plt.cla()
    plt.imshow(ising_array)
    plt.scatter(bug_coord[1], bug_coord[0], marker = 'o', s = 0.5, color = 'yellow')
    plt.pause(0.00001)


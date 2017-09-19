import sys, os
import numpy as np
import random
from math import exp
import matplotlib.pyplot as plt
import time
import cProfile
import re




n=200
ising_array = [[-1 if (random.random() > 0.5) else 1  for i in range(n)] for j in range(n)]
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
            energy += -J*ising_array[i][(j+1)%n] - J*ising_array[(i+1)%n][j] + H*ising_array[i][j]
    return energy


def calc_del_e(spin, J):
    global ising_array
    i = spin[0]
    j = spin[1]
    neigh = [ising_array[(i+1)%n][j], ising_array[i][(j+1)%n], ising_array[(i-1)%n][j], ising_array[i][(j-1)%n]]
    #when ising array is a numpy array neigh = ising_array[[(i+1)%n,i,(i-1)%n,i],[j,(j+1)%n,j,(j-1)%n]]
    sigma_spin = ising_array[i][j]
    old_e = -J*sum([sigma_spin*i for i in neigh]) + H*sigma_spin
    new_e = -J*sum([-sigma_spin*i for i in neigh]) - H*sigma_spin
    return (new_e - old_e)


print(calc_total_energy(J))

tot_en = calc_total_energy(J)

random.seed(1002)
sweeps = 500

for s in range(sweeps):
    for sw in range(n*n):
        spin = [random.randint(0,n-1),random.randint(0,n-1)]
        del_e = calc_del_e(spin, J)
        if(exp(-beta*del_e) > random.random()):
            ising_array[spin[0]][spin[1]] *= -1
            tot_en += del_e
    plt.cla()
    plt.imshow(ising_array)
    plt.pause(0.0001)


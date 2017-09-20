import sys, os
import numpy as np
import random
from math import exp
import matplotlib.pyplot as plt
import time
import cProfile
import re



random.seed(13829)

n=200
n_bugs = 50
ising_array = [[-1 if (random.random() > 0.5) else 1  for i in range(n)] for j in range(n)]
bug_coord = [[random.randint(0,n-1) for i in range(n_bugs)] for j in range(2)]
J = 0
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
n_bug_steps = 20
alive = [1 for i in range(n_bugs)]
#initial probabilities for moving forwards, backwards or turning
pf = 0.25
pb = 0.25
pt = 0.5
###
#pt = 1
#pf = 0
#pb = 0
prev_magnetization = 0
shift_J = 0
for s in range(sweeps):
    #update ising base
    prev_J = J
    J = J + 0.1*(random.random()-0.5 + shift_J)
    for sw in range(n*n):
        spin = [random.randint(0,n-1),random.randint(0,n-1)]
        del_e = calc_del_e(spin, J)
        if(exp(-beta*del_e) > random.random()):
            ising_array[spin[0]][spin[1]] *= -1
            tot_en += del_e
    #move bugs
    food_consumed = [0 for i in range(n_bugs)]
    prev_step = [[0 if (random.random() < 0.5) else 1 for i in range(n_bugs)]]
    prev_step.append([1 if prev_step[0][i] == 0 else 0 for i in range(n_bugs)])
    move_stats = np.zeros((n_bugs,3))
    for bsw in range(n_bug_steps):
        prev = list(bug_coord)
        for cbug in range(n_bugs):
            #move random step forward backward or turning
            p0 = bug_coord[0][cbug]
            p1 = bug_coord[1][cbug]
            choose_dir = random.random()
            if choose_dir < pf:
                #move forwards
                bug_coord[0][cbug] += prev_step[0][cbug]
                bug_coord[1][cbug] += prev_step[1][cbug]
                move_stats[cbug][0] += 1
            elif choose_dir < pf + pb:
                #move backwards
                bug_coord[0][cbug] -= prev_step[0][cbug]
                bug_coord[1][cbug] -= prev_step[1][cbug]
                move_stats[cbug][1] += 1
            else:
                #turn
                move_stats[cbug][2] += 1
                if random.random() < 0.5:
                    #turn left
                    bug_coord[0][cbug] -= prev_step[1][cbug]
                    bug_coord[1][cbug] += prev_step[0][cbug]
                else:
                    #turn right
                    bug_coord[0][cbug] += prev_step[1][cbug]
                    bug_coord[1][cbug] -= prev_step[0][cbug]
            prev_step[0][cbug] = bug_coord[0][cbug] - p0#prev[0][cbug]
            prev_step[1][cbug] = bug_coord[1][cbug] - p1#prev[1][cbug]
            bug_coord[0][cbug] = bug_coord[0][cbug]%n
            bug_coord[1][cbug] = bug_coord[1][cbug]%n
            #eat up!
            if ising_array[bug_coord[0][cbug]][bug_coord[1][cbug]] == 1:
                food_consumed[cbug] += 1
                ising_array[bug_coord[0][cbug]][bug_coord[1][cbug]] = -1
    mean_food = np.mean(np.asarray(food_consumed))
    #is mean best?
    if mean_food > 0:
        successful_stats = np.ones(3)
        for ci,i in enumerate(food_consumed):
            if i > mean_food:
                successful_stats += move_stats[ci]
        norm = np.sum(successful_stats)
        pf = successful_stats[0]/norm
        pb = successful_stats[1]/norm
        pt = successful_stats[2]/norm
    #change so that all moves have some probability
    #track_move_probs.append([pf, pb, pt])
    #
    #
    #
    # change ising params:
    magnetization = 0
    for i in ising_array:
        for j in i:
            magnetization += j
    print(magnetization - prev_magnetization)
    if magnetization - prev_magnetization > 0:
        shift_J = 0.05*(J - prev_J) #reinforce change
    else:
        shift_J = 0.05*(prev_J - J)
    prev_magnetization = magnetization
    print(pf, pb, pt, J, magnetization)
    plt.cla()
    plt.imshow(ising_array)
    plt.scatter(bug_coord[1], bug_coord[0], marker = 'o', s = 5, color = 'red')
    plt.pause(0.00001)


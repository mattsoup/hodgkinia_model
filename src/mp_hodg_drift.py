#!/usr/bin/env python


import multiprocessing as mp
import random
import sys
import time

import numpy as np
import scipy.stats

num_insects = 100        #insect population size
num_hodg = 10            #hodgkinia bottleneck size
adult_hodg = 3000        #hodgkinia population size in adult insects
mut_rate = 100            #Inverse of the mutation rate
num_generations = 2001    #Number of generations
inflection_point = 0.7    #Inflection point on the sigmoidal fitness curve
num_genes = 10            #Number of genes per Hodgkinia genome

#First populate the insect population, and each insect with a hodgkinia population,
#and each hodgkinia with functional genes
insect_pop = []
for x in range(0, num_insects):
    insect_pop.append([])
    for y in range(0, num_hodg): # The number of appends equals the number of genes
        insect_pop[x].append([])
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)
        insect_pop[x][y].append(1)

#Populates a list of host fitnesses, which for now are equal
fitness_list = []
for x in range(0, num_insects):
    fitness_list.append(1 / num_insects)

#Function that 'grows' the Hodgkinia population from the bottleneck
#size to the adult size
def hodg_growth(insect_pop):
    for x in range(0, num_insects):
        for y in range(0, num_hodg):
            first = insect_pop[x][y][0]
            second = insect_pop[x][y][1]
            third = insect_pop[x][y][2]
            fourth = insect_pop[x][y][3]
            fifth = insect_pop[x][y][4]
            sixth = insect_pop[x][y][5]
            seventh = insect_pop[x][y][6]
            eighth = insect_pop[x][y][7]
            ninth = insect_pop[x][y][8]
            tenth = insect_pop[x][y][9]
            for z in range(1, (adult_hodg // num_hodg)): # starts at one because there's already one there
                insect_pop[x].append([first, second, third, fourth, fifth, sixth, seventh, eighth, ninth, tenth])
    return insect_pop

#Function that calls 10 processes t run the 'single mutation' function
def mp_mutation(insect_pop):
    processes = []
    prev_end = 0
    temp = []
    for x in range(10):
        start_index = prev_end
        end_index = ((num_insects // 10) * (x + 1))
        processes.append(mp.Process(target=single_mutation,
                                    args=(start_index, end_index, temp)))
        prev_end = ((num_insects // 10) * (x + 1))
    for p in processes:
        p.daemon = True
        p.start()
    for p in processes:
        while p.is_alive() == True:
            time.sleep(0.1)
        p.join()
    insect_pop = []
    for x in range(0, 100, 10):
        temp_file = open("%s_%s.out" % (x, (x + 10)), "r")
        temp = []
        new_insect = []
        for line in temp_file:
            if line.startswith("new"):
                if len(new_insect) != 0:
                    insect_pop.append(new_insect)
                new_insect = []
                temp = []
            else:
                for gene in line[:-1]:
                    temp.append(int(gene))
                new_insect.append(temp)
                temp = []
        temp_file.close()
    return insect_pop

#Function that mutates the Hodgkinia genes
def single_mutation(start_index, end_index, temp):
    temp_out = open("%s_%s.out" % (start_index, end_index), "w")
    slice = insect_pop[start_index:end_index]
    for x in range(0, len(slice)):
        for y in range(0, len(slice[x])):
            for z in range(0, num_genes):
                mut = random.randint(0, mut_rate)
                if mut == 7:
                    #print "MUTANT", x, y, z
                    if slice[x][y][z] == 1:
                        slice[x][y][z] = 0
                    else:
                        pass
    for item in slice:
        temp_out.write("new\n")
        for hodg in item:
            for gene in hodg:
                temp_out.write("%s" % gene)
            temp_out.write("\n")
    temp_out.write("new\n")
    temp_out.close()
    return


#Function to reproduce the insect population, based on their fitnesses
def insect_reproduction(insect_pop, fitness_list):
    new_insect_pop = []
    fitness_list = fitness_list
    #Make a list of indices from the insect population, randomly chosen based
    #on weight (fitness)
    one_temp = np.random.choice(len(insect_pop), num_insects, p=fitness_list)
    #Populate a temporary list with what will be the new insect population
    temp = []
    for item in one_temp:
        temp.append(insect_pop[item])
    #Populate the new insect population
    for x in range(0, num_insects):
        new_insect_pop.append([])
        for y in range(0, num_hodg):
            temp_hodg = random.choice(temp[y])
            first = temp_hodg[0]
            second = temp_hodg[1]
            third = temp_hodg[2]
            fourth = temp_hodg[3]
            fifth = temp_hodg[4]
            sixth = temp_hodg[5]
            seventh = temp_hodg[6]
            eighth = temp_hodg[7]
            ninth = temp_hodg[8]
            tenth = temp_hodg[9]
            new_insect_pop[x].append([first, second, third, fourth, fifth, sixth, seventh, eighth, ninth, tenth])
    #Calculate the proportion of each gene that has been lost in each insect
    for x in range(0, num_insects):
        gene1 = 0
        gene2 = 0
        gene3 = 0
        gene4 = 0
        gene5 = 0
        gene6 = 0
        gene7 = 0
        gene8 = 0
        gene9 = 0
        gene10 = 0
        lost_gene1 = 0
        lost_gene2 = 0
        lost_gene3 = 0
        lost_gene4 = 0
        lost_gene5 = 0
        lost_gene6 = 0
        lost_gene7 = 0
        lost_gene8 = 0
        lost_gene9 = 0
        lost_gene10 = 0
        for y in range(0, num_hodg):
            gene1 += 1
            gene2 += 1
            gene3 += 1
            gene4 += 1
            gene5 += 1
            gene6 += 1
            gene7 += 1
            gene8 += 1
            gene9 += 1
            gene10 += 1
            if new_insect_pop[x][y][0] == 0:
                lost_gene1 += 1
            if new_insect_pop[x][y][1] == 0:
                lost_gene2 += 1
            if new_insect_pop[x][y][2] == 0:
                lost_gene3 += 1
            if new_insect_pop[x][y][3] == 0:
                lost_gene4 += 1
            if new_insect_pop[x][y][4] == 0:
                lost_gene5 += 1
            if new_insect_pop[x][y][5] == 0:
                lost_gene6 += 1
            if new_insect_pop[x][y][6] == 0:
                lost_gene7 += 1
            if new_insect_pop[x][y][7] == 0:
                lost_gene8 += 1
            if new_insect_pop[x][y][8] == 0:
                lost_gene9 += 1
            if new_insect_pop[x][y][9] == 0:
                lost_gene10 += 1
        #"proportion_geneX" is the proportion of that gene that has been lost
        proportion_gene1 = lost_gene1 / gene1
        proportion_gene2 = lost_gene2 / gene2
        proportion_gene3 = lost_gene3 / gene3
        proportion_gene4 = lost_gene4 / gene4
        proportion_gene5 = lost_gene5 / gene5
        proportion_gene6 = lost_gene6 / gene6
        proportion_gene7 = lost_gene7 / gene7
        proportion_gene8 = lost_gene8 / gene8
        proportion_gene9 = lost_gene9 / gene9
        proportion_gene10 = lost_gene10 / gene10
        #Multiply together the proportions of each gene that has been lost
        #times =  (1 - proportion_gene1) * (1 - proportion_gene2) * (1 - proportion_gene3) * (1 - proportion_gene4) * (1 - proportion_gene5) * (1 - proportion_gene6) * (1 - proportion_gene7) * (1 - proportion_gene7) * (1 - proportion_gene8) * (1 - proportion_gene10)
        #Find the average proportion of each gene lost
        #avg = ((1 - proportion_gene1) + (1 - proportion_gene2) + (1 - proportion_gene3) + (1 - proportion_gene4) + (1 - proportion_gene5) + (1 - proportion_gene6) + (1 - proportion_gene7) + (1 - proportion_gene8) + (1 - proportion_gene9) + (1 - proportion_gene10)) / float(10)
        #Find the harmonic mean of the proportion of genes retained
        hmean_list = [1-proportion_gene1, 1-proportion_gene2, 1-proportion_gene3, 1-proportion_gene4, 1-proportion_gene5, 1-proportion_gene6, 1-proportion_gene7, 1-proportion_gene8, 1-proportion_gene9, 1-proportion_gene10]
        #Find the fitness for the insect
        #Force the fitness to be zero if all copies of any gene has been lost
        if 0 in hmean_list:
            fitness = 0
        else:
            hmean = scipy.stats.hmean(hmean_list)
            fitness = (1 / (1+(np.exp(-2*((hmean * num_hodg) - (inflection_point * num_hodg))))))
        fitness_list[x] = fitness

    fitness_sum = sum(fitness_list)
    avg_fitness = (fitness_sum / len(fitness_list))
    fitness_range = max(fitness_list) - min(fitness_list)
    for y in range(0, len(fitness_list)):
        fitness_list[y] /= fitness_sum
    return new_insect_pop, fitness_list, avg_fitness, fitness_range

out = open(sys.argv[1], "w")
out.write("Mutation rate: %s\n" % mut_rate)
out.write("Generations: %s\n" % num_generations)
out.write("Inflection point: %s\n\n" % inflection_point)
out_genotypes = open(sys.argv[1] + ".genotypes", "w")

out.write("Generation\tTotal genes\tLost genes\tAverage fitness\tRange of fitnesses\tCooperators\tNine\tEight\tSeven\tSix\tFive\tFour\tThree\tTwo\tOne\tSelfish\n")
#Runs the model for X number of generations, keeps track of genes lost, etc.
for x in range(num_generations):
    print("Generation %s" % (x + 1))
    insect_pop = hodg_growth(insect_pop)
    insect_pop = mp_mutation(insect_pop)
    insect_pop, fitness_list, avg_fitness, fitness_range = insect_reproduction(insect_pop, fitness_list)
    lost_genes = 0
    total_genes = 0
    cooperators = 0
    fragmented = 0
    selfish = 0
    nine = 0
    eight = 0
    seven = 0
    six = 0
    five = 0
    four = 0
    three = 0
    two = 0
    one = 0
    for y in range(0, len(insect_pop)):
        for z in range(0, len(insect_pop[y])):
            total_genes += 10
            my_sum = sum(insect_pop[y][z])
            if my_sum == 10:
                cooperators += 1
            elif my_sum == 0:
                selfish += 1
            if my_sum == 9:
                nine += 1
            if my_sum == 8:
                eight += 1
            if my_sum == 7:
                seven += 1
            if my_sum == 6:
                six += 1
            if my_sum == 5:
                five += 1
            if my_sum == 4:
                four += 1
            if my_sum == 3:
                three += 1
            if my_sum == 2:
                two += 1
            if my_sum == 1:
                one += 1
            else:
                lost_genes += (10 - my_sum)
    print("Total genes: %s" % total_genes)
    print("Lost_genes: %s" % lost_genes)
    print("Average fitness: %s" % avg_fitness)
    print("Range of fitnesses: %s" % fitness_range)
    print("Cooperators: %s" % cooperators)
    print("Nine: %s" % nine)
    print("Eight: %s" % eight)
    print("Seven: %s" % seven)
    print("Six: %s" % six)
    print("Five: %s" % five)
    print("Four: %s" % four)
    print("Three: %s" % three)
    print("Two: %s" % two)
    print("One: %s" % one)
    print("Selfish: %s\n" % selfish)
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (x + 1, total_genes, lost_genes, avg_fitness, fitness_range, cooperators, nine, eight, seven, six, five, four, three, two, one, selfish))
    out_genotypes.write("Generation %s\n" % (x + 1))
    for y in range(len(insect_pop)):
        out_genotypes.write("Insect %s:\n" % (y + 1))
        for symbiont in insect_pop[y]:
            out_genotypes.write("%s\n" % symbiont)
    if lost_genes == total_genes:
        out.write("All genes lost after %s generations" % (x + 1))
        print("All genes lost after %s generations" % (x + 1))
        break

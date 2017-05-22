#!/usr/bin/env python
'''Simulation of Hodgkinia and symbionts
'''

import multiprocessing as mp
import os
import random
import sys
import time

import numpy as np
import scipy.stats

num_insects = 100        #insect population size
num_hodg = 11            #hodgkinia bottleneck size
adult_hodg_factor = 300  #factor by which to increase the num_hodg by in adult insects
adult_hodg = num_hodg * adult_hodg_factor        #hodgkinia population size in adult insects
mut_rate = 100           #Inverse of the mutation rate
num_generations = 2001   #Number of generations
inflection_point = 0.7   #Inflection point on the sigmoidal fitness curve
num_genes = 10           #Number of genes per Hodgkinia genome
mutation_mean = (num_insects * adult_hodg * num_genes) / mut_rate  #Average number of mutations to introduce each insect generation

#some stupid variables for the new mutation function
cicada_len = len(str(num_insects)) if len(str(num_insects)) == len(str(num_insects - 1)) else len(str(num_insects - 1))
hodg_len = len(str(adult_hodg)) if len(str(adult_hodg)) == len(str(adult_hodg - 1)) else len(str(adult_hodg - 1))
gene_len = len(str(num_genes)) if len(str(num_genes)) == len(str(num_genes - 1)) else len(str(num_genes - 1))

output_dir = sys.argv[1]
try:
    os.mkdir(output_dir)
except FileExistsError:
    if not os.path.isdir(output_dir):
        print("A file with that name already exists")
        sys.exit(-1)

#First populate the insect population, and each insect with a hodgkinia population,
#and each hodgkinia with functional genes
insect_pop = []
for insect in range(0, num_insects):
    insect_pop.append([])
    for hodg in range(0, num_hodg):
        insect_pop[insect].append([])
        insect_pop[insect][hodg].extend([1] * num_genes)

#Populates a list of host fitnesses, which for now are equal
fitness_list = []
for insect in range(0, num_insects):
    fitness_list.append(1 / float(num_insects)) # one of these needs to be a float, else they will all be zero


def hodg_growth(my_insect_pop):
    '''Function that 'grows' the Hodgkinia population from the bottleneck
       size to the adult size'''
    for x in range(0, num_insects):
        for y in range(0, num_hodg):
            my_genes = []
            for my_gene in range(num_genes):
                my_genes.append(my_insect_pop[x][y][my_gene])
            for _z in range(1, (adult_hodg // num_hodg)):
                # starts at one because there's already one there
                my_insect_pop[x].append(list(my_genes))
    return my_insect_pop

def mp_mutation(my_insect_pop):
    '''Function that calls 10 processes t run the 'single mutation' function'''
    processes = []
    prev_end = 0
    temp = []
    for x in range(10):
        start_index = prev_end
        end_index = ((num_insects // 10) * (x + 1))
        processes.append(mp.Process(target=single_mutation,
                                    args=(start_index, end_index)))
        prev_end = ((num_insects // 10) * (x + 1))
    for p in processes:
        p.daemon = True
        p.start()
    for p in processes:
        while p.is_alive():
            time.sleep(0.1)
        p.join()
    my_insect_pop = []
    for x in range(0, 100, 10):
        temp_file = open(output_dir + os.sep + "%s_%s.out" % (x, (x + 10)), "r")
        temp = []
        new_insect = []
        for line in temp_file:
            if line.startswith("new"):
                if len(new_insect) != 0:
                    my_insect_pop.append(new_insect)
                new_insect = []
                temp = []
            else:
                for my_gene in line[:-1]:
                    temp.append(int(my_gene))
                new_insect.append(temp)
                temp = []
        temp_file.close()
    return my_insect_pop


def single_mutation(start_index, end_index):
    '''Function that mutates the Hodgkinia genes'''
    temp_out = open(output_dir + os.sep + "%s_%s.out" % (start_index, end_index), "w")
    my_slice = insect_pop[start_index:end_index]
    for x in range(0, len(my_slice)):
        for y in range(0, len(my_slice[x])):
            for z in range(0, num_genes):
                mut = random.randint(0, mut_rate)
                if mut == 7:
                    #print "MUTANT", x, y, z
                    if my_slice[x][y][z] == 1:
                        my_slice[x][y][z] = 0
                    else:
                        pass
    for item in my_slice:
        temp_out.write("new\n")
        for my_hodg in item:
            for my_gene in my_hodg:
                temp_out.write("%s" % my_gene)
            temp_out.write("\n")
    temp_out.write("new\n")
    temp_out.close()
    return

def all_mutations(insect_pop):
	'''Function that mutates Hodgkinia genes much faster.
	   Currently will only work if hodg_len < cicada_len <= gene_len'''
	num_mutations = int(np.random.normal(mutation_mean, mutation_mean / 4))
	mutants = random.sample(range(0, mutation_mean * mut_rate), num_mutations)
	for mutation in mutants:
		mutation = str(mutation)
		while len(mutation) < len(str(mutation_mean * mut_rate)):
			mutation = "0" + mutation
		hodg = mutation[0:hodg_len]
		cicada = mutation[hodg_len:hodg_len + cicada_len]
		gene = mutation[-gene_len:]
		if insect_pop[int(cicada)][int(hodg)][int(gene)] == 1:
			insect_pop[int(cicada)][int(hodg)][int(gene)] = 0
	return insect_pop
		
def insect_reproduction(my_insect_pop, my_fitness_list):
    '''Function to reproduce the insect population, based on their fitnesses'''
    new_insect_pop = []
    #Make a list of indices from the insect population, randomly chosen based
    #on weight (fitness)
    one_temp = np.random.choice(len(my_insect_pop),
                                num_insects, p=my_fitness_list)
    #Populate a temporary list with what will be the new insect population
    temp = []
    for item in one_temp:
        temp.append(my_insect_pop[item])
    #Populate the new insect population
    for x in range(0, num_insects):
        new_insect_pop.append([])
        for y in range(0, num_hodg):
            temp_hodg = random.choice(temp[y])
            new_insect_pop[x].append(temp_hodg)
    #Calculate the proportion of each gene that has been lost in each insect
    for x in range(0, num_insects):
        genes = [0] * num_genes
        my_lost_genes = [0] * num_genes
        for y in range(0, num_hodg):
            genes = list(map(lambda x: x + 1, range(num_genes)))
            for my_gene in range(num_genes):
                my_lost_genes[my_gene] += (1 - new_insect_pop[x][y][0])
        #"proportion_geneX" is the proportion of that gene that has been lost
        proportion_genes = []
        for my_gene in range(num_genes):
            proportion_genes.append(my_lost_genes[my_gene] / genes[my_gene])
        #Multiply together the proportions of each gene that has been lost
        #times =  (1 - proportion_gene1) * (1 - proportion_gene2) * (1 - proportion_gene3) * (1 - proportion_gene4) * (1 - proportion_gene5) * (1 - proportion_gene6) * (1 - proportion_gene7) * (1 - proportion_gene7) * (1 - proportion_gene8) * (1 - proportion_gene10)
        #Find the average proportion of each gene lost
        #avg = ((1 - proportion_gene1) + (1 - proportion_gene2) + (1 - proportion_gene3) + (1 - proportion_gene4) + (1 - proportion_gene5) + (1 - proportion_gene6) + (1 - proportion_gene7) + (1 - proportion_gene8) + (1 - proportion_gene9) + (1 - proportion_gene10)) / float(10)
        #Find the harmonic mean of the proportion of genes retained
        hmean_list = list(map(lambda x: 1 - x, proportion_genes))
        #Find the fitness for the insect
        #Force the fitness to be zero if all copies of any gene has been lost
        if 0 in hmean_list:
            fitness = 0
        else:
            hmean = scipy.stats.hmean(hmean_list)
            fitness = (1 / (1+(np.exp(-2*((hmean * num_hodg) - (inflection_point * num_hodg))))))
        my_fitness_list[x] = fitness

    fitness_sum = sum(my_fitness_list)
    my_avg_fitness = (fitness_sum / len(my_fitness_list))
    my_fitness_range = max(my_fitness_list) - min(my_fitness_list)
    for y in range(0, len(my_fitness_list)):
        my_fitness_list[y] /= fitness_sum
    return new_insect_pop, my_fitness_list, my_avg_fitness, my_fitness_range

out = open(output_dir + os.sep + "results.txt", "w")
out.write("Mutation rate: %s\n" % mut_rate)
out.write("Generations: %s\n" % num_generations)
out.write("Inflection point: %s\n\n" % inflection_point)
out_genotypes = open(output_dir + os.sep + "results.genotypes", "w")

out.write("\t".join(["Generation", "Total genes", "Lost genes",
                     "Average fitness", "Range of fitnesses",
                     "Cooperators", "Nine", "Eight", "Seven", "Six",
                     "Five", "Four", "Three", "Two", "One", "Selfish"]))
out.write("\n")
#Runs the model for X number of generations, keeps track of genes lost, etc.
for generation in range(num_generations):
    start = time.time()
    print("Generation %s" % (generation + 1))
    insect_pop = hodg_growth(insect_pop)
    insect_pop = all_mutations(insect_pop)
    #insect_pop = mp_mutation(insect_pop)
    insect_pop, fitness_list, avg_fitness, fitness_range = \
        insect_reproduction(insect_pop, fitness_list)
    total_genes = 0
    fragmented = 0
    activated_genes = [0] * (num_genes + 1)
    lost_genes = 0
    for insect in range(0, len(insect_pop)):
        for gene in range(0, len(insect_pop[insect])):
            total_genes += num_genes
            my_sum = sum(insect_pop[insect][gene])
            activated_genes[my_sum] += 1
            lost_genes += num_genes - my_sum
    print("Total genes: %s" % total_genes)
    print("Lost_genes: %s" % lost_genes)
    print("Average fitness: %s" % avg_fitness)
    print("Range of fitnesses: %s" % fitness_range)
    print("Cooperators: %s" % activated_genes[num_genes])
    for cnt_gene in range(num_genes, 0, -1):
        print('%d:\t%d' % (cnt_gene, activated_genes[cnt_gene]))
    print("Selfish: %s\n" % activated_genes[0])
    out.write("\t".join(map(lambda x: str(x),
                            [generation + 1, total_genes, lost_genes,
                             avg_fitness, fitness_range] +
                            activated_genes)))
    out.write('\n')
    out_genotypes.write("Generation %s\n" % (generation + 1))
    for num_insect, insect in enumerate(insect_pop):
        out_genotypes.write("Insect %s:\n" % (num_insect + 1))
        for symbiont in insect:
            out_genotypes.write("%s\n" % symbiont)
    if lost_genes == total_genes:
        out.write("All genes lost after %s generations" % (generation + 1))
        print("All genes lost after %s generations" % (generation + 1))
        break
    end = time.time()
    print end - start

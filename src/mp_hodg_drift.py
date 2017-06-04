#!/usr/bin/env python
'''Simulation of Hodgkinia and symbionts
'''

import json
import multiprocessing as mp
import os
import random
import sys
import time

import numpy as np
import scipy.stats

conf = json.load(open(sys.argv[1]))

# This is a bad piece of code, do not emulate
num_genes = mut_rate = num_insects = adult_hodg_factor = num_hodg = None
num_generations = None
for key, val in conf.items():
    op = "float" if val.find('.') > -1 else "int"
    print("%s = %s(%s)" % (key, op, val))
    exec("%s = %s(%s)" % (key, op, val))

output_dir = sys.argv[2]
start = time.time()

# hodgkinia population size in adult insects
adult_hodg = num_hodg * adult_hodg_factor

# Average number of mutations to introduce each insect generation
mutation_mean = (num_insects * adult_hodg * num_genes) / float(mut_rate)


try:
    os.mkdir(output_dir)
except FileExistsError:
    if not os.path.isdir(output_dir):
        print("A file with that name already exists")
        sys.exit(-1)

#First populate the insect population, and each insect with a hodgkinia population,
#and each hodgkinia with functional genes
insect_pop = [1] * (num_insects * num_hodg * num_genes)
insect_pop[7] = 0
#Populates a list of host fitnesses, which for now are equal
fitness_list = [(1 / float(num_insects))] * num_insects     # one of these needs to be a float, else they will all be zero

class hodg_lineage():
    parent = []
    genotype = []
    origin_generation = ""
    origin_hodg = ""

lineage_list = []
for x in range(num_insects * num_hodg):
    my_lineage = hodg_lineage()
    my_lineage.parent = "1" * num_genes
    my_lineage.genotype = "1" * num_genes
    my_lineage.origin_generation = "0"
    my_lineage.origin_hodg = "0"
    lineage_list.append(my_lineage)

def hodg_growth(my_insect_pop, my_lineage_list):
    '''Function that 'grows' the Hodgkinia population from the bottleneck
       size to the adult size'''
    new_insect_pop = []
    new_lineage_list = []
    for x in range(len(my_lineage_list)):
    	for y in range(adult_hodg_factor):
    		new_lineage_list.append(my_lineage_list[x])
    for x in range(0, len(my_insect_pop), num_genes):
        for y in range(adult_hodg_factor):
        	for z in range(num_genes):
	        	new_insect_pop.append(my_insect_pop[x + z])
    return new_insect_pop, new_lineage_list

def all_mutations(my_insect_pop, my_lineage_list):
    '''Function that mutates Hodgkinia genes much faster.'''
    num_mutations = int(np.random.normal(mutation_mean, mutation_mean / 20)) # Picks a number of mutations, normally distributed around mutation_mean
	#mutants = random.sample(range(0, mutation_mean * mut_rate), num_mutations)
    mutants = [random.randint(0,(num_insects * adult_hodg * num_genes) - 1) for x in range(num_mutations)]
    for mutant in mutants:
        if my_insect_pop[mutant] == 1:
            current_hodg = mutant / num_genes
            parent = "".join(str(x) for x in my_insect_pop[current_hodg * num_genes:(current_hodg * num_genes) + num_genes])
            my_insect_pop[mutant] = 0
            genotype = "".join(str(x) for x in my_insect_pop[current_hodg * num_genes:(current_hodg * num_genes) + num_genes])
            my_lineage = hodg_lineage()
            my_lineage.parent = parent + "," + str(my_lineage_list[current_hodg].origin_generation) + "," + str(my_lineage_list[current_hodg].origin_hodg)
            my_lineage.genotype = genotype
            my_lineage.origin_generation = generation + 1
            my_lineage.origin_hodg = current_hodg
            my_lineage_list[current_hodg] = my_lineage
    return my_insect_pop, my_lineage_list


def insect_reproduction(my_insect_pop, my_fitness_list, my_lineage_list):
    '''Function to reproduce the insect population, based on their fitnesses'''
    new_insect_pop = []
    new_lineage_list = []
    #Make a list of indices from the insect population, randomly chosen based
    #on weight (fitness)
    to_reproduce = np.random.choice(num_insects,
                                num_insects, p=my_fitness_list)
    #Populate a temporary list with what will be the new insect population
    #Populate the new insect population
    for item in to_reproduce:
        my_start = item * adult_hodg * num_genes
        my_end = my_start + (adult_hodg * num_genes)
        for x in range(num_hodg):
	        my_hodg = random.randint(0, adult_hodg - 1)
	        new_lineage_list.append(my_lineage_list[(item * adult_hodg) + my_hodg])
	        for y in range(my_start + (my_hodg * num_genes), my_start + (my_hodg * num_genes) + num_genes):
	            new_insect_pop.append(my_insect_pop[y])
    #Calculate the proportion of each gene that has been lost in each insect
    for x in range(num_insects):
        my_start = x * num_hodg * num_genes
        my_end = my_start + (num_hodg * num_genes)
        genes = [0] * num_genes
        my_lost_genes = [0] * num_genes
        for y in range(my_start, my_end):
            my_gene = (y + num_genes) % num_genes
            if new_insect_pop[y] == 0:
                my_lost_genes[my_gene] += 1
#         for y in range(0, num_hodg):
#             genes = list(map(lambda x: x + 1, range(num_genes)))
#             for z in range(num_genes):
#                 my_lost_genes[z] += (1 - new_insect_pop[x][y][z])
        #"proportion_geneX" is the proportion of that gene that has been lost
        proportion_genes = []
        for my_gene in range(num_genes):
            proportion_genes.append(my_lost_genes[my_gene] / float(num_hodg)) #One needs to be a float else there will be many ones and zeroes
        #Multiply together the proportions of each gene that has been lost
        #times =  (1 - proportion_geinsect_pne1) * (1 - proportion_gene2) * (1 - proportion_gene3) * (1 - proportion_gene4) * (1 - proportion_gene5) * (1 - proportion_gene6) * (1 - proportion_gene7) * (1 - proportion_gene7) * (1 - proportion_gene8) * (1 - proportion_gene10)
        #Find the harmonic mean of the proportion of genes retained
        hmean_list = list(map(lambda x: 1 - x, proportion_genes))
        #print hmean_list
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
    return new_insect_pop, my_fitness_list, my_avg_fitness, my_fitness_range, new_lineage_list

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
    print("Generation %s" % (generation + 1))
    insect_pop, lineage_list = hodg_growth(insect_pop, lineage_list)
    insect_pop, lineage_list = all_mutations(insect_pop, lineage_list)
    insect_pop, fitness_list, avg_fitness, fitness_range, lineage_list = \
        insect_reproduction(insect_pop, fitness_list, lineage_list)
    total_genes = 0
    fragmented = 0
    active_genes = [0] * (num_genes + 1)
    lost_genes = 0
    for insect in range(num_insects):
        my_start = insect * num_hodg * num_genes
        my_end = my_start + (num_hodg * num_genes)
        for my_hodg in range(my_start, my_end, num_genes):
            my_genes = insect_pop[my_hodg:my_hodg + num_genes]
            my_sum = sum(my_genes)
            active_genes[num_genes - my_sum] += 1
            total_genes += num_genes
            lost_genes += (num_genes - my_sum)
    print("Total genes: %s" % total_genes)
    print("Lost_genes: %s" % lost_genes)
    print("Average fitness: %s" % avg_fitness)
    print("Range of fitnesses: %s" % fitness_range)
    #print("Cooperators: %s" % active_genes[num_genes])
    for cnt_gene in range(num_genes + 1):
        print('%d genes:\t%d' % (num_genes - cnt_gene, active_genes[cnt_gene]))
    print("\n")
    out.write("\t".join(map(lambda x: str(x),
                            [generation + 1, total_genes, lost_genes,
                             avg_fitness, fitness_range] +
                            active_genes)))
    out.write('\n')
    out_genotypes.write("Generation %s\n" % (generation + 1))
    for x in range(len(lineage_list)):
        for item in lineage_list[x].parent:
        	out_genotypes.write(str(item))
        out_genotypes.write("\t")
        for item in lineage_list[x].genotype:
        	out_genotypes.write(str(item))
        out_genotypes.write("\t")
        out_genotypes.write("%s\t%s\n" % (lineage_list[x].origin_generation, lineage_list[x].origin_hodg))
        
    if lost_genes == total_genes:
        out.write("All genes lost after %s generations" % (generation + 1))
        print("All genes lost after %s generations" % (generation + 1))
        break
end = time.time()
print("This script took %s seconds" % (end - start))

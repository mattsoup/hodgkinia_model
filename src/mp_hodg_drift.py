#!/usr/bin/env python
'''Simulation of Hodgkinia and symbionts
'''

import copy
import json
import multiprocessing as mp
import os
import random
import sys
import time

import numpy as np
import scipy.stats


class hodg_lineage():
    parent = []
    genotype = []
    origin_generation = ""
    origin_hodg = ""


def hodg_growth(my_insect_pop, my_lineage_list, adult_hodg_factor,
                num_hodg, num_insects, num_genes):
    '''Function that 'grows' the Hodgkinia population from the bottleneck
       size to the adult size'''
    new_insect_pop = []
    new_lineage_list = []
#     print(my_insect_pop)
    #for x in range(len(my_lineage_list)):
    #    for y in range(adult_hodg_factor):
    #        new_lineage_list.append(my_lineage_list[x])

    # Grows all Hodgkinia cells by a set number (adult_hodg_factor)
#     for factor in range(adult_hodg_factor):
#         new_lineage_list.extend(copy.copy(my_lineage_list))
#         new_insect_pop.extend(copy.copy(my_insect_pop))
    to_reproduce = [np.random.choice(num_hodg, (adult_hodg_factor * num_hodg))
                    for x in range(num_insects)]
#     print(to_reproduce)
    for x in range(0, len(to_reproduce)):
        insect_start = x * num_hodg * num_genes
        insect_end = insect_start + (num_hodg * num_genes)
#         print("insect", x)
        for my_hodg in to_reproduce[x]:
            genes = my_insect_pop[insect_start + (my_hodg * num_genes):insect_start + (my_hodg * num_genes) + num_genes]
#             print(genes, my_hodg)
            for gene in genes:
                new_insect_pop.append(gene)
#     print new_insect_pop
    #for x in range(0, len(my_insect_pop), num_genes):
    #    for y in range(adult_hodg_factor):
    #        for z in range(num_genes):
    #            new_insect_pop.append(my_insect_pop[x + z])
    return new_insect_pop, new_lineage_list


def all_mutations(my_insect_pop, my_lineage_list, mutation_mean,
                  num_insects, adult_hodg, num_genes,
                  generation):
    '''Function that mutates Hodgkinia genes much faster.'''
    num_mutations = int(np.random.normal(mutation_mean, mutation_mean / 20))
    # Picks a number of mutations, normally distributed around mutation_mean
    mutants = [random.randint(0, (num_insects * adult_hodg * num_genes) - 1)
               for x in range(num_mutations)]
    for mutant in mutants:
        if my_insect_pop[mutant] == 1:
            current_hodg = mutant // num_genes
#             parent = "".join(
#                 str(x)
#                 for x in my_insect_pop[current_hodg * num_genes:
#                                        (current_hodg * num_genes) + num_genes])
            my_insect_pop[mutant] = 0
#             genotype = "".join(
#                 str(x)
#                 for x in my_insect_pop[current_hodg * num_genes:
#                                        (current_hodg * num_genes) + num_genes])
#             my_lineage = hodg_lineage()
#             my_lineage.parent = parent + "," + str(
#                 my_lineage_list[current_hodg].origin_generation) + "," + str(
#                     my_lineage_list[current_hodg].origin_hodg)
#             my_lineage.genotype = genotype
#             my_lineage.origin_generation = generation + 1
#             my_lineage.origin_hodg = current_hodg
#             my_lineage_list[current_hodg] = my_lineage
    return my_insect_pop, my_lineage_list


def insect_reproduction(my_insect_pop, my_fitness_list, my_lineage_list,
                        num_insects, adult_hodg, num_genes, num_hodg,
                        inflection_point, k, low_fitness):
    '''Function to reproduce the insect population, based on their fitnesses'''
    new_insect_pop = []
    new_lineage_list = []
    # Make a list of indices from the insect population, randomly chosen based
    # on weight (fitness)
    my_sum = float(sum(my_fitness_list))
    for x in range(len(my_fitness_list)):
        my_fitness_list[x] = my_fitness_list[x] / my_sum
    to_reproduce = np.random.choice(num_insects,
                                    num_insects, p=my_fitness_list)
    # Populate a temporary list with what will be the new insect population
    # Populate the new insect population
    for item in to_reproduce:
        my_start = item * adult_hodg * num_genes
        my_end = my_start + (adult_hodg * num_genes)
        for x in range(num_hodg):
            my_hodg = random.randint(0, adult_hodg - 1)
#             new_lineage_list.append(my_lineage_list[(item * adult_hodg) +
#                                                     my_hodg])
            for y in range(my_start + (my_hodg * num_genes), my_start +
                           (my_hodg * num_genes) + num_genes):
                new_insect_pop.append(my_insect_pop[y])
    # Calculate the proportion of each gene that has been lost in each insect
    for x in range(num_insects):
        my_start = x * num_hodg * num_genes
        my_end = my_start + (num_hodg * num_genes)
        genes = [0] * num_genes
        my_lost_genes = [0] * num_genes
        for y in range(my_start, my_end):
            my_gene = (y + num_genes) % num_genes
            if new_insect_pop[y] == 0:
                my_lost_genes[my_gene] += 1
        proportion_genes = []
        for my_gene in range(num_genes):
            proportion_genes.append(my_lost_genes[my_gene] / float(num_hodg))
            # One needs to be a float else there will be many ones and zeroes
        # Find the harmonic mean of the proportion of genes retained
        hmean_list = list(map(lambda x: 1 - x, proportion_genes))
        my_fitness_list[x] = calc_fitness(hmean_list, k, inflection_point, low_fitness)

    fitness_sum = sum(my_fitness_list)
    #assert (fitness_sum != 0), "All insects are extinct"
    my_avg_fitness = fitness_sum / len(my_fitness_list)
    my_fitness_range = max(my_fitness_list) - min(my_fitness_list)
    for y in range(0, len(my_fitness_list)):
        my_fitness_list[y] /= fitness_sum
    return (new_insect_pop, my_fitness_list, my_avg_fitness,
            my_fitness_range, new_lineage_list)

def calc_fitness(my_hmean_list, k, inflection_point, low_fitness):
    if 0 not in my_hmean_list:
        hmean = scipy.stats.hmean(my_hmean_list)
        fitness = (1 / (1 + (np.exp(-k * (hmean - 
                   inflection_point)))))
    else:
        fitness = low_fitness
    return fitness

class Conf:
    def __init__(self, conf_file):
        my_conf = json.load(open(conf_file))
        for key, val in my_conf.items():
            op = float if val.find('.') > -1 else int
            print("%s = %s(%s)" % (key, op, val))
            setattr(self, key, op(val))


def simulate(conf_file, output_dir, silent=True):
    c = Conf(conf_file)

    # hodgkinia population size in adult insects
    adult_hodg = c.num_hodg * c.adult_hodg_factor

    # Average number of mutations to introduce each insect generation
    mutation_mean = float(c.num_insects * adult_hodg * c.num_genes) / c.mut_rate
    
    low_fitness = c.low_fitness

    try:
        os.mkdir(output_dir)
    except FileExistsError:
        if not os.path.isdir(output_dir):
            print("A file with that name already exists")
            return None

    lineage_list = []
    for x in range(c.num_insects * c.num_hodg):
        my_lineage = hodg_lineage()
        my_lineage.parent = "1" * c.num_genes
        my_lineage.genotype = "1" * c.num_genes
        my_lineage.origin_generation = "0"
        my_lineage.origin_hodg = "0"
        lineage_list.append(my_lineage)

    # First populate the insect population, and each insect with a hodgkinia
    # population, and each hodgkinia with functional genes
    insect_pop = [1] * (c.num_insects * c.num_hodg * c.num_genes)
    # Populates a list of host fitnesses, which for now are equal
    fitness_list = [float(1) / c.num_insects] * c.num_insects
    # one of these needs to be a float, else they will all be zero
    out = open(output_dir + os.sep + "results.txt", "w")
    out.write("Mutation rate: %s\n" % c.mut_rate)
    out.write("Generations: %s\n" % c.num_generations)
    out.write("Inflection point: %s\n\n" % c.inflection_point)
    out_genotypes = open(output_dir + os.sep + "results.genotypes", "w")

    out.write("\t".join(["Generation", "Total genes", "Lost genes",
                         "Average fitness", "Range of fitnesses",
                         "Cooperators", "Nine", "Eight", "Seven", "Six",
                         "Five", "Four", "Three", "Two", "One", "Selfish"]))
    out.write("\n")
    # Runs the model for X number of generations
    # keeps track of genes lost, etc.
    for generation in range(c.num_generations):
#         if not silent:
#             print("Generation %s" % (generation + 1))
        insect_pop, lineage_list = hodg_growth(insect_pop,
                                               lineage_list,
                                               c.adult_hodg_factor,
                                               c.num_hodg,
                                               c.num_insects,
                                               c.num_genes)
        insect_pop, lineage_list = all_mutations(insect_pop,
                                                 lineage_list, mutation_mean,
                                                 c.num_insects,
                                                 adult_hodg,
                                                 c.num_genes,
                                                 generation)
        insect_pop, fitness_list, avg_fitness, fitness_range, lineage_list = \
            insect_reproduction(insect_pop, fitness_list, lineage_list,
                                c.num_insects, adult_hodg, c.num_genes,
                                c.num_hodg, c.inflection_point, c.k, low_fitness)
        total_genes = 0
        fragmented = 0
        active_genes = [0] * (c.num_genes + 1)
        lost_genes = 0
        for insect in range(c.num_insects):
            my_start = insect * c.num_hodg * c.num_genes
            my_end = my_start + c.num_hodg * c.num_genes
            for my_hodg in range(my_start, my_end, c.num_genes):
                my_genes = insect_pop[my_hodg:my_hodg + c.num_genes]
                my_sum = sum(my_genes)
                active_genes[c.num_genes - my_sum] += 1
                total_genes += c.num_genes
                lost_genes += (c.num_genes - my_sum)
#         if not silent:
#             print("Total genes: %s" % total_genes)
#             print("Lost_genes: %s" % lost_genes)
#             print("Average fitness: %s" % avg_fitness)
#             print("Range of fitnesses: %s" % fitness_range)
#             # print("Cooperators: %s" % active_genes[num_genes])
#             for cnt_gene in range(c.num_genes + 1):
#                 print('%d genes:\t%d' % (c.num_genes - cnt_gene,
#                       active_genes[cnt_gene]))
#             print("\n")

        out.write("\t".join(map(lambda x: str(x),
                                [generation + 1, total_genes, lost_genes,
                                 avg_fitness, fitness_range] +
                                active_genes)))
        out.write('\n')
#         out_genotypes.write("Generation %s\n" % (generation + 1))
#         for x in range(len(lineage_list)):
#             for item in lineage_list[x].parent:
#                 out_genotypes.write(str(item))
#             out_genotypes.write("\t")
#             for item in lineage_list[x].genotype:
#                 out_genotypes.write(str(item))
#             out_genotypes.write("\t")
#             out_genotypes.write("%s\t%s\n" % (
#                 lineage_list[x].origin_generation,
#                 lineage_list[x].origin_hodg))

        if lost_genes == total_genes:
            out.write("All genes lost after %s generations" % (generation + 1))
            if not silent:
                print("All genes lost after %s generations" % (generation + 1))
            break

if __name__ == "__main__":
    start = time.time()
    conf_file = sys.argv[1]
    out_dir = sys.argv[2]
    simulate(conf_file, out_dir, False)
    end = time.time()
    print("This script took %s seconds" % (end - start))

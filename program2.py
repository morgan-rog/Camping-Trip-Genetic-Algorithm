from random import *
import numpy as np
import time
import copy

class Thing:
    def __init__(self, utility, weight):
        self.utility = utility
        self.weight = weight

    def __str__(self):
        return f'Utility: {self.utility} Weight: {self.weight}'


class Genome:
    def __init__(self, binary_string):
        self.binary_string = binary_string
        self.fitness_value = 0

    def set_fitness(self, g_fitness):
        self.fitness_value = g_fitness
    
    def __str__(self):
        return f'genome: {self.binary_string} - fitness: {self.fitness_value}'


def generate_genome(length):
    binary_string = choices([0, 1], k=length)
    new_genome = Genome(binary_string)
    return new_genome


def generate_population(size, genome_length):
    return [generate_genome(genome_length) for i in range(size)]


def fitness(genome, my_things, weight_limit):
    if len(genome.binary_string) != len(my_things):
        raise ValueError('genome and my_things must be same length')

    weight = 0
    value = 0

    for i, thing in enumerate(my_things):
        if genome.binary_string[i] == 1:
            weight += thing.weight
            value += thing.utility

        if weight > weight_limit:
            return 1

    return value


def L2(population):
    fitness_list = []

    for g in population:
        fitness_list.append(g.fitness_value)

    fitness_array = np.array(fitness_list)
    fitness_norm = np.linalg.norm(fitness_array)
    fitness_normalized = fitness_array / fitness_norm
    return fitness_normalized


def selection_pair(population, L2_norm):
    a_genome, b_genome = choices(population, weights=L2_norm, k=2)
    return a_genome, b_genome


def single_point_crossover(a_genome, b_genome):
    '''cuts two genomes at a random point and swaps them to return a tuple of two new genomes'''

    if len(a_genome.binary_string) != len(b_genome.binary_string):
        raise ValueError('genomes a and b must be the same length')

    length = len(a_genome.binary_string)
    if length < 2:
        return a_genome.binary_string, b_genome.binary_string

    # random point to cut between 1 and the length of the genome minus 1
    p = randint(1, length - 1)
    a_genome.binary_string = a_genome.binary_string[0:p] + b_genome.binary_string[p:]
    b_genome.binary_string = b_genome.binary_string[0:p] + a_genome.binary_string[p:]

    return a_genome, b_genome


def mutation(genome, probability = 0.0001):
    for index in range(len(genome.binary_string)):
        
        # random() returns random float uniformly in semi-open range [0.0, 1.0]
        if random() < probability:
            # genome bit gets flipped
            # abs(0 - 1) = 1 and abs(1 - 1) = 0
            genome.binary_string[index] = abs(genome.binary_string[index] - 1)

    return genome


def calculate_average_fitness(population):
    total_fitness = 0
    for genome in population:
        total_fitness += genome.fitness_value
    avg_fitness = total_fitness / len(population)
    return avg_fitness


def get_things_from_file(input_file):
    things_list = []

    with open(input_file, 'r') as file:
        for line in file.readlines():
            utility, weight = line.split("\t")

            utility = float(utility.strip())
            weight = float(weight.strip())
            thing = Thing(utility, weight)
            things_list.append(thing)

    return things_list


def choose_random_things(num, the_things):
    things = the_things[:]
    random_things = []
    for i in range(num):
        random_thing = choices(things, weights=[item.utility for item in things], k=1)[0]
        random_things.append(random_thing)
        things.remove(random_thing)

    return random_things


def record_avg_fitness(avg_fitness):
    with open('average_fitness.txt', 'a') as outfile:
        outfile.write(str(avg_fitness))
        outfile.write('\n')


def record_max_fitness_genome(max_fitness_genome, my_things):
    things_taken = []
    total_utility = 0
    total_weight = 0

    for i, thing in enumerate(my_things):
        if max_fitness_genome.binary_string[i] == 1:
            things_taken.append(thing)
            total_utility += thing.utility
            total_weight += thing.weight

    with open('max_fitness_genome.txt', 'a') as outfile:
        outfile.write('max fitness genome:\n\t')
        outfile.write(str(max_fitness_genome))
        outfile.write(f'\n\ttotal weight: {total_weight}')
        outfile.write(f'\n\ttotal utility: {total_utility}')
        outfile.write('\n\nItems taken:\n\t')
        for thing in things_taken:
            outfile.write(str(thing))
            outfile.write("\n\t")
        outfile.write('\n\n')

            
def clear_fitness_file():
    f = open('average_fitness.txt', 'w')
    f.close()


def clear_fitness_genome_file():
    f = open('max_fitness_genome.txt', 'w')
    f.close()


if __name__ == '__main__':
    weight_limit = 500
    genome_and_things_size = 20 # size of genome (bits) and amount of items
    population_size_list = [500, 1000, 1500, 2000, 2500]
    number_of_generations = 1000
    iteration = 50
    all_things = get_things_from_file('Program2Input.txt')
    # clear output files
    clear_fitness_file() 
    clear_fitness_genome_file()

    for pop_size in population_size_list:
        my_things = choose_random_things(genome_and_things_size, all_things)
        population = generate_population(pop_size, genome_and_things_size)
        max_fitness = 0
        max_fitness_genome = copy.copy(population[0])

        start_time = time.time() 

        # calculate fitness for each genome in the population
        for g in population:
            g.set_fitness(fitness(g, my_things, weight_limit))

        total_gen_count = 1
        gen_count_iteration = 1
        last_iter_avg_fitness = calculate_average_fitness(population)

        # output average fitness to average_fitness.txt
        record_avg_fitness(last_iter_avg_fitness)

        for i in range(number_of_generations - 1):
            # sort population based on fitness function for each genome
            population = sorted(population, key=lambda genome: genome.fitness_value, reverse=True)

            generation_max_fitness_genome = copy.copy(population[0])
            generation_max_fitness = population[0].fitness_value
            if generation_max_fitness > max_fitness:
                max_fitness = generation_max_fitness
                max_fitness_genome = generation_max_fitness_genome

            # next generation starts with two genomes from population that have the highest fitness scores
            next_generation = population[0:2]

            # use L2 normalization to convert fitness scores to probablility distribution
            L2_normalized = L2(population)
            l = np.linalg.norm(L2_normalized)

            # create offspring and add it to next generation
            for j in range(int(len(population) / 2) - 1):
                parents = selection_pair(population, L2_normalized)
                offspring_a, offspring_b = single_point_crossover(parents[0], parents[1])
                offspring_a = mutation(offspring_a)
                offspring_a.set_fitness(fitness(offspring_a, my_things, weight_limit))
                offspring_b = mutation(offspring_b)
                offspring_b.set_fitness(fitness(offspring_b, my_things, weight_limit))
                next_generation += [offspring_a, offspring_b]

            population = next_generation

            if gen_count_iteration == iteration:
                gen_count_iteration = 0
                one_percent_last_avg = last_iter_avg_fitness * .01
                new_avg_fitness = calculate_average_fitness(population)
                if (abs(new_avg_fitness - last_iter_avg_fitness)) < one_percent_last_avg:
                    break
                else:
                    avg_fitness = new_avg_fitness
                    # record new last iteration average fitness for next iteration
                    last_iter_avg_fitness = new_avg_fitness
            else:
                avg_fitness = calculate_average_fitness(population)

            # output average fitness to average_fitness.txt
            record_avg_fitness(avg_fitness)

            total_gen_count += 1
            gen_count_iteration += 1
        

        end_time = time.time()
        total_time = round(end_time - start_time, 2)

        record_max_fitness_genome(max_fitness_genome, my_things)

        print(f'Starting population: {pop_size}')
        print(f'Max fitness after {total_gen_count} generations: {round(max_fitness, 2)}')
        print(f'Average fitness: {round(calculate_average_fitness(population), 2)}')
        print(f'time: {total_time}s')
        print()
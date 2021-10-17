# given input text file containing 400 sets of (utility, weight) pairs
# utility: float ranged 0-10. 0 or 1 is something we could do without, 9-10 is vital
# task is to use genetic algorithm to find good selection of items to pack while staying within weight guidelines (< 500 lbs)

# use initial population of 1000 random selections
# use mutation rate of 0.0001 -> each item in each selection has, independently, a 1/10,000 chance of being changed during mutation step
# use L2 normalization to convert fitness scores to probability distribution
# assign fitness score of 1 when it exceeds 500 lb limit
# record average fitness for each generation and save to an output file. Continue iterating until average fitness improves less than 1% across 10 generations
# report highest fitness selection found: items taken and what total utility is -> put this in a text file

# prepare short report (1 or 2 pages) describing program -> data structures used and why, any places able to parallelize or optimize program, overview of what running the program was like
# (overall efficiency, how many generations needed, etc.)
from random import *

class Thing:
    def __init__(self, utility, weight):
        self.utility = utility
        self.weight = weight

    def __str__(self):
        return f'Utility: {self.utility} Weight: {self.weight}'


def generate_genome(length):
    return choices([0, 1], k=length)


def generate_population(size, genome_length):
    return [generate_genome(genome_length) for i in range(size)]


def fitness(genome, my_things, weight_limit):
    # if len(genome) != len(my_things):
    #     raise ValueError("genome and my_things must be same length")

    weight = 0
    value = 0

    for i, thing in enumerate(my_things):
        if genome[i] == 1:
            weight += thing.weight
            value += thing.utility

        if weight > weight_limit:
            return 0

    return value


def selection_pair(population, my_things, weight_limit):
    a_genome, b_genome = choices(population, weights=[fitness(genome, my_things, weight_limit) for genome in population], k=2)
    return a_genome, b_genome


def single_point_crossover(a_genome, b_genome):
    '''cuts two genomes at a random point and swaps them to return a tuple of two new genomes'''

    if len(a_genome) != len(b_genome):
        raise ValueError("genomes a and b must be the same length")

    length = len(a_genome)
    if length < 2:
        return a_genome, b_genome

    # random point to cut between 1 and the length of the genome minus 1
    p = randint(1, length - 1)
    return a_genome[0:p] + b_genome[p:], b_genome[0:p] + a_genome[p:]


def mutation(genome, probability = 0.0001):
    for index in range(len(genome)):
        
        # random() returns random float uniformly in semi-open range [0.0, 1.0]
        if random() < probability:
            # abs(0 - 1) = 1
            # abs(1 - 1) = 0
            # genome bit gets flipped
            genome[index] = abs(genome[index] - 1)


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


def choose_random_things(num, things):
    random_things = []
    for i in range(num):
        thing = choice(things)
        random_things.append(thing)

    return random_things



# MAIN
weight_limit = 500
genome_and_things_size = 20 # size of genome (in bits) and amount of items
population_size = 500
number_of_generations = 200
things = get_things_from_file('Program2Input.txt')

my_things = choose_random_things(genome_and_things_size, things)
population = generate_population(population_size, genome_and_things_size)
max_fitness_list = []

for i in range(number_of_generations):
    # sort population based on fitness function for each genome
    print(my_things)
    population = sorted(population, key=lambda genome: fitness(genome, my_things, weight_limit), reverse=True)
    max_fitness = fitness(population[0], my_things, weight_limit)
    max_fitness_list.append(max_fitness)
    next_generation = population[0:2]

    for j in range(int(len(population) / 2) - 1):
        parents = selection_pair(population, my_things, weight_limit)
        offspring_a, offspring_b = single_point_crossover(parents[0], parents[1])
        offspring_a = mutation(offspring_a)
        offspring_b = mutation(offspring_b)
        next_generation += [offspring_a, offspring_b]

    population = next_generation

sorted_max_fitness = sorted(max_fitness_list, reverse=True)
max_fitness = sorted_max_fitness[0]
print(f'max fitness for {number_of_generations} generations: {max_fitness}')
    



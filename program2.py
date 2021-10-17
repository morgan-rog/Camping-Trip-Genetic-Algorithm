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


def selection_pair(population):
    a_genome, b_genome = choices(population, weights=[genome.fitness_value for genome in population], k=2)
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
    print(f'total fitness: {total_fitness}')
    print(f'population size: {len(population)}')
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


# MAIN
weight_limit = 500
genome_and_things_size = 5 # size of genome (bits) and amount of items
population_size = 1000
number_of_generations = 100
all_things = get_things_from_file('Program2Input.txt')

my_things = choose_random_things(genome_and_things_size, all_things)
population = generate_population(population_size, genome_and_things_size)
max_fitness = 0

# calculate fitness for each genome in the population
for g in population:
    g.set_fitness(fitness(g, my_things, weight_limit))


for i in range(number_of_generations):
    # sort population based on fitness function for each genome
    population = sorted(population, key=lambda genome: genome.fitness_value, reverse=True)

    generation_max_fitness = population[0].fitness_value
    if generation_max_fitness > max_fitness:
        max_fitness = generation_max_fitness

    # next generation starts with two genomes from population that have the highest fitness scores
    next_generation = population[0:2]

    # create offspring and add it to next generation
    for j in range(int(len(population) / 2) - 1):
        parents = selection_pair(population)
        offspring_a, offspring_b = single_point_crossover(parents[0], parents[1])
        offspring_a = mutation(offspring_a)
        offspring_a.set_fitness(fitness(offspring_a, my_things, weight_limit))
        offspring_b = mutation(offspring_b)
        offspring_b.set_fitness(fitness(offspring_b, my_things, weight_limit))
        next_generation += [offspring_a, offspring_b]


    population = next_generation

print(f'max fitness after {number_of_generations} generations: {max_fitness}')
print(f'average fitness: {calculate_average_fitness(population)}')
    



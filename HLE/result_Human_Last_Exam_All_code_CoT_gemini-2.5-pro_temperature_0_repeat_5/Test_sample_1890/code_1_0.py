import random

def simulate_evolution(population_size, gene_length, generations, mutation_rate, strain_name):
    """
    Simulates the emergence of mutations in a population.

    - Frequency is the count of individuals with any mutation.
    - Spectrum is the count of unique mutation patterns.
    """
    # Start with a population of identical, non-mutated individuals
    # A gene is represented by a list of 0s. A mutation turns a 0 to a 1.
    population = [[0] * gene_length for _ in range(population_size)]
    original_gene = tuple([0] * gene_length) # Use tuple for set hashing

    for _ in range(generations):
        for i in range(population_size):
            # Give each gene a chance to mutate
            for j in range(gene_length):
                if random.random() < mutation_rate:
                    population[i][j] = 1 # A mutation occurs at this position

    # After all generations, analyze the final population
    mutated_individuals = 0
    mutated_gene_types = set()

    for gene in population:
        gene_tuple = tuple(gene)
        if gene_tuple != original_gene:
            mutated_individuals += 1
            mutated_gene_types.add(gene_tuple)

    # Print the results for this strain
    print(f"--- {strain_name} Strain Simulation ---")
    print(f"Parameters: Population={population_size}, Generations={generations}, Mutation Rate={mutation_rate}")
    print(f"Result: Increased Frequency of Variants?")
    print(f"A total of {mutated_individuals} out of {population_size} individuals became variants.")
    print(f"Result: Wider Spectrum of Mutations?")
    print(f"A total of {len(mutated_gene_types)} unique mutation patterns (spectrum) were found.")
    print("-" * 35 + "\n")

# --- Simulation Parameters ---
# A hypermutator can have a 100x to 1000x higher mutation rate.
POPULATION = 1000
GENE_LENGTH = 50 # Represents the mucA gene
GENERATIONS = 10
WILD_TYPE_RATE = 0.0001
HYPERMUTATOR_RATE = 0.01 # 100x higher than wild-type

# --- Run Simulations ---
# This simulation demonstrates the principle described in the biological context.
# A higher mutation rate leads to more variants and a wider variety of mutations.
simulate_evolution(POPULATION, GENE_LENGTH, GENERATIONS, WILD_TYPE_RATE, "Wild-Type")
simulate_evolution(POPULATION, GENE_LENGTH, GENERATIONS, HYPERMUTATOR_RATE, "Hypermutator")
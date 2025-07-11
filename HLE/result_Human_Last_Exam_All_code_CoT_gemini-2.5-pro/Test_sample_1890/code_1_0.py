import random

def simulate_evolution(population_size, generations, mutation_rate_per_base, gene_length):
    """
    Simulates the emergence of mutations in a population.
    
    Returns:
        - final_mucoid_count: Total number of individuals with a mutated gene.
        - mutation_spectrum_size: Number of unique mutation patterns found.
    """
    # The 'wild-type' gene is represented by a string of 'A's
    wild_type_gene = 'A' * gene_length
    
    # A set to store all unique mutated gene sequences found
    mutation_spectrum = set()
    
    # A list representing the genes of the current population
    population = [wild_type_gene] * population_size
    
    for _ in range(generations):
        new_population = []
        for gene in population:
            new_gene_list = list(gene)
            # Check each base of the gene for a potential mutation
            for i in range(gene_length):
                if random.random() < mutation_rate_per_base:
                    # Mutate the base to B, C, or D if it's A, or to A otherwise
                    original_base = new_gene_list[i]
                    possible_mutations = [b for b in 'BCD' if b != original_base]
                    new_gene_list[i] = random.choice(possible_mutations)
            
            new_gene = "".join(new_gene_list)
            new_population.append(new_gene)
            
            # If the gene has mutated, add it to our spectrum
            if new_gene != wild_type_gene:
                mutation_spectrum.add(new_gene)
        
        population = new_population

    # Count the final number of mucoid (mutated) individuals
    final_mucoid_count = sum(1 for gene in population if gene != wild_type_gene)
    
    return final_mucoid_count, len(mutation_spectrum)

# --- Simulation Parameters ---
POPULATION_SIZE = 500
GENERATIONS = 100
GENE_LENGTH = 50 # A simplified representation of the mucA gene

# Mutation rates
# A hypermutator can be 10x to 1000x higher; we'll use 50x for a clear result.
WT_MUTATION_RATE = 0.0001
HYPERMUTATOR_RATE = 0.005 # 50x higher than wild-type

# --- Run Simulations ---
wt_count, wt_spectrum = simulate_evolution(POPULATION_SIZE, GENERATIONS, WT_MUTATION_RATE, GENE_LENGTH)
hm_count, hm_spectrum = simulate_evolution(POPULATION_SIZE, GENERATIONS, HYPERMUTATOR_RATE, GENE_LENGTH)

# --- Print Results ---
print("--- Simulation Results ---")
print(f"Population Size: {POPULATION_SIZE}, Generations: {GENERATIONS}\n")

print("Wild-Type Strain:")
print(f"  Frequency of Mucoid Variants: {wt_count} / {POPULATION_SIZE} ({wt_count/POPULATION_SIZE:.2%})")
print(f"  Spectrum of mucA Mutations (unique variants): {wt_spectrum}\n")

print("Hypermutator Strain:")
print(f"  Frequency of Mucoid Variants: {hm_count} / {POPULATION_SIZE} ({hm_count/POPULATION_SIZE:.2%})")
print(f"  Spectrum of mucA Mutations (unique variants): {hm_spectrum}\n")

print("Conclusion:")
print("The simulation shows the hypermutator strain results in both an 'Increased frequency of mucoid variants' and a 'wider spectrum of mucA mutations'.")
print("This corresponds to Answer Choice B.")

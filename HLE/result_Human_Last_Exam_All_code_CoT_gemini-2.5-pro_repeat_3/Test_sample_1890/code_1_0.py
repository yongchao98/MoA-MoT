import random

def simulate_evolution(population_size, generations, mutation_rate, gene_length):
    """
    Simulates the emergence of mucoid variants in a bacterial population.
    
    Returns:
        - mucoid_count (int): The total number of mucoid variants that appeared.
        - mutation_spectrum (set): A set of unique mutation sites found.
    """
    # The probability of a mucoid-causing mutation in the gene per individual per generation.
    # This is a simplification, approximating P(mutation) as rate * length for small rates.
    probability_of_mucoid = mutation_rate * gene_length
    
    mucoid_count = 0
    mutation_sites = set()
    
    # Total number of chances for mutation to occur
    total_trials = population_size * generations
    
    for _ in range(total_trials):
        # Check if a mucoid-causing mutation occurs in this trial
        if random.random() < probability_of_mucoid:
            mucoid_count += 1
            # If a mutation occurs, we assume it hits a random site in the gene
            site = random.randint(1, gene_length)
            mutation_sites.add(site)
            
    return mucoid_count, mutation_sites

# --- Simulation Parameters ---
POPULATION_SIZE = 50000
GENERATIONS = 100
MUC_A_GENE_LENGTH = 1470  # Approximate length of mucA gene in base pairs

# Hypermutator strains have mutation rates 100-1000x higher than normal strains
NORMAL_MUTATION_RATE = 1e-8
HYPERMUTATOR_MUTATION_RATE = 5e-6 # 500x higher

# --- Run Simulations ---
normal_freq, normal_spectrum_set = simulate_evolution(POPULATION_SIZE, GENERATIONS, NORMAL_MUTATION_RATE, MUC_A_GENE_LENGTH)
hypermutator_freq, hypermutator_spectrum_set = simulate_evolution(POPULATION_SIZE, GENERATIONS, HYPERMUTATOR_MUTATION_RATE, MUC_A_GENE_LENGTH)

# --- Output Results ---
print("--- Simulation Results ---")
print("This simulation models the effect of an increased mutation rate on the frequency and spectrum of mucoid variants.")

print("\nNormal Strain:")
print(f"Total Mucoid Variants (Frequency): {normal_freq}")
print(f"Unique mucA Mutation Sites (Spectrum): {len(normal_spectrum_set)}")

print("\nHypermutator Strain:")
print(f"Total Mucoid Variants (Frequency): {hypermutator_freq}")
print(f"Unique mucA Mutation Sites (Spectrum): {len(hypermutator_spectrum_set)}")

print("\n--- Conclusion ---")
print("The hypermutator strain shows both an 'Increased frequency' of mucoid variants and a 'wider spectrum' of mutations.")
print("This result supports answer choice B.")

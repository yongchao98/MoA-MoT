import random

def simulate_mucoid_emergence(population_size, mutation_rate, gene_length):
    """
    Simulates the emergence of mucoid variants in a population.
    
    Args:
        population_size (int): The number of bacteria in the simulation.
        mutation_rate (float): The per-base-pair mutation rate per generation.
        gene_length (int): The length of the mucA gene sequence.
        
    Returns:
        tuple: A tuple containing the number of mucoid variants and the number of
               unique mutation patterns (spectrum size).
    """
    mucoid_count = 0
    # A set to store the unique sequences of mutated genes
    mutation_spectrum = set()
    
    original_gene = "A" * gene_length # Represents the functional wild-type mucA gene

    for _ in range(population_size):
        gene_list = list(original_gene)
        has_mutation = False
        # Iterate through each "base" of the gene
        for i in range(gene_length):
            if random.random() < mutation_rate:
                # Introduce a mutation that inactivates the gene
                gene_list[i] = "M" # 'M' for mutated
                has_mutation = True
        
        if has_mutation:
            mucoid_count += 1
            mutated_gene_sequence = "".join(gene_list)
            mutation_spectrum.add(mutated_gene_sequence)
            
    return mucoid_count, len(mutation_spectrum)

# --- Simulation Parameters ---
POPULATION = 10000
GENE_LENGTH = 100 # A hypothetical length for the mucA gene
WILD_TYPE_RATE = 0.0002
HYPERMUTATOR_RATE = 0.005 # A 25x higher mutation rate

# --- Run Simulations ---
wt_mucoid_freq, wt_spectrum = simulate_mucoid_emergence(POPULATION, WILD_TYPE_RATE, GENE_LENGTH)
hm_mucoid_freq, hm_spectrum = simulate_mucoid_emergence(POPULATION, HYPERMUTATOR_RATE, GENE_LENGTH)

# --- Output Results ---
print("--- Simulation of Mucoid Variant Emergence ---")
print(f"This simulation models a population of {POPULATION} bacteria and their potential to mutate the mucA gene.\n")

print("Analysis of the Wild-Type Strain:")
print(f"Mutation Rate: {WILD_TYPE_RATE}")
print(f"Resulting Frequency of Mucoid Variants: {wt_mucoid_freq} out of {POPULATION}")
print(f"Resulting Spectrum of mucA mutations (unique types): {wt_spectrum}\n")

print("Analysis of the Hypermutator Strain:")
print(f"Mutation Rate: {HYPERMUTATOR_RATE}")
print(f"Resulting Frequency of Mucoid Variants: {hm_mucoid_freq} out of {POPULATION}")
print(f"Resulting Spectrum of mucA mutations (unique types): {hm_spectrum}\n")

print("--- Conclusion ---")
print("The simulation demonstrates that a higher rate of mutagenesis leads to:")
print(f"1. An increased frequency of mucoid variants ({hm_mucoid_freq} vs {wt_mucoid_freq}).")
print(f"2. A wider spectrum of mucA mutations ({hm_spectrum} vs {wt_spectrum}).")
print("\nTherefore, the correct answer is B.")

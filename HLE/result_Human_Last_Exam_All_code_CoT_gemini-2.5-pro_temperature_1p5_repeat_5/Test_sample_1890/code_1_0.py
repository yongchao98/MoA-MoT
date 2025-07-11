import random

def simulate_mutations(population_size, gene_length, mutation_rate_per_base):
    """
    Simulates mutations in a gene for a bacterial population.

    Returns a tuple containing:
    - The number of mucoid variants (frequency).
    - The number of unique mucoid genotypes (spectrum).
    """
    # The 'original' gene is represented by a string of 'A's.
    original_gene = ['A'] * gene_length
    
    mucoid_genotypes = []

    for _ in range(population_size):
        new_gene = list(original_gene)
        is_mutated = False
        for i in range(gene_length):
            # Check if a random mutation occurs at this base.
            if random.random() < mutation_rate_per_base:
                new_gene[i] = 'M' # 'M' for mutated base
                is_mutated = True
        
        # If the gene was mutated, this bacterium is a mucoid variant.
        if is_mutated:
            mucoid_genotypes.append("".join(new_gene))
            
    # Frequency is the total count of mucoid variants.
    mucoid_frequency = len(mucoid_genotypes)
    # Spectrum is the count of unique mutated gene sequences.
    mucoid_spectrum = len(set(mucoid_genotypes))
    
    return mucoid_frequency, mucoid_spectrum

def main():
    """Main function to run the simulation and print results."""
    # Simulation parameters
    POPULATION_SIZE = 10000
    GENE_LENGTH = 1000 # Length of the mucA gene
    NORMAL_MUTATION_RATE = 0.00001 # A low, baseline mutation rate
    HYPERMUTATOR_RATE = 0.0005 # A 50x higher mutation rate

    print("Running simulation to model the effect of hypermutation...")
    print(f"Population Size: {POPULATION_SIZE}")
    print(f"Normal Mutation Rate: {NORMAL_MUTATION_RATE}")
    print(f"Hypermutator Mutation Rate: {HYPERMUTATOR_RATE}\n")

    # Run simulation for the normal population
    normal_freq, normal_spectrum = simulate_mutations(POPULATION_SIZE, GENE_LENGTH, NORMAL_MUTATION_RATE)

    # Run simulation for the hypermutator population
    hypermutator_freq, hypermutator_spectrum = simulate_mutations(POPULATION_SIZE, GENE_LENGTH, HYPERMUTATOR_RATE)

    print("--- Simulation Results ---")
    print("Normal Population:")
    print(f"  Frequency of mucoid variants: {normal_freq} out of {POPULATION_SIZE}")
    print(f"  Spectrum of mucA mutations: {normal_spectrum} unique variants")
    
    print("\nHypermutator Population:")
    print(f"  Frequency of mucoid variants: {hypermutator_freq} out of {POPULATION_SIZE}")
    print(f"  Spectrum of mucA mutations: {hypermutator_spectrum} unique variants")

    print("\n--- Conclusion ---")
    print("The simulation shows that the hypermutator strain's increased rate of mutagenesis leads to:")
    print(f"1. An 'Increased frequency' of mucoid variants ({hypermutator_freq} is much greater than {normal_freq}).")
    print(f"2. A 'wider spectrum' of mucA mutations ({hypermutator_spectrum} is much greater than {normal_spectrum}).")
    print("\nThis directly corresponds to Answer Choice B.")


if __name__ == "__main__":
    main()

<<<B>>>
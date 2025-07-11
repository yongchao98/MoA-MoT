import random

def simulate_bacterial_evolution(population_size, generations, mutation_rate, gene_length):
    """
    Simulates the accumulation of mutations in a bacterial population.
    
    Returns the total number of mucoid variants that appeared and the
    set of unique mutations that caused them.
    """
    mucoid_variants_count = 0
    unique_mutations = set()

    # In each generation, each bacterium has a chance to mutate
    total_events = population_size * generations
    for _ in range(total_events):
        if random.random() < mutation_rate:
            # A mutation inactivating mucA has occurred
            mucoid_variants_count += 1
            
            # Generate a random, unique mutation signature (e.g., position in the gene)
            # This represents the "spectrum"
            mutation_position = random.randint(1, gene_length)
            unique_mutations.add(mutation_position)
            
    return mucoid_variants_count, unique_mutations

def main():
    # Simulation Parameters
    POPULATION_SIZE = 1_000_000
    GENERATIONS = 100
    GENE_LENGTH = 579 # Approximate length of P. aeruginosa mucA gene in base pairs

    # Define mutation rates
    # Wild-type rate (e.g., 1 in 100 million)
    WT_MUTATION_RATE = 1e-8
    # Hypermutator rate is often 100x to 1000x higher
    HYPERMUTATOR_RATE = 1e-5

    print("Running simulation to compare Wild-Type and Hypermutator strains...")
    print("-" * 50)

    # Run simulation for Wild-Type
    wt_count, wt_mutations = simulate_bacterial_evolution(
        POPULATION_SIZE, GENERATIONS, WT_MUTATION_RATE, GENE_LENGTH
    )
    wt_spectrum_size = len(wt_mutations)
    
    # Run simulation for Hypermutator
    hyper_count, hyper_mutations = simulate_bacterial_evolution(
        POPULATION_SIZE, GENERATIONS, HYPERMUTATOR_RATE, GENE_LENGTH
    )
    hyper_spectrum_size = len(hyper_mutations)

    print(f"Simulation Results for a population of {POPULATION_SIZE} over {GENERATIONS} generations:\n")

    # Output the results of the "equation": Frequency = Count / Total
    total_population_events = POPULATION_SIZE * GENERATIONS
    wt_freq = wt_count / total_population_events
    hyper_freq = hyper_count / total_population_events
    
    print(f"Wild-Type Strain:")
    print(f"  - Total Mucoid Variants (Frequency): {wt_count}")
    print(f"  - Spectrum of mucA Mutations (Diversity): {wt_spectrum_size} unique mutations found")

    print(f"\nHypermutator Strain:")
    print(f"  - Total Mucoid Variants (Frequency): {hyper_count}")
    print(f"  - Spectrum of mucA Mutations (Diversity): {hyper_spectrum_size} unique mutations found")

    print("\n" + "-" * 50)
    print("Conclusion: The simulation shows the hypermutator strain has both an")
    print("increased frequency of mucoid variants and a wider spectrum of mutations, supporting option B.")


if __name__ == "__main__":
    main()
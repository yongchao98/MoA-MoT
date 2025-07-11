import random

def run_bacterial_evolution_simulation():
    """
    Simulates the emergence of mucoid variants in wild-type and hypermutator
    populations of Pseudomonas aeruginosa.
    """
    # --- Parameters ---
    # The mucA gene is represented as a string. A mutation changes a character.
    GENE_LENGTH = 30
    WILD_TYPE_GENE = 'A' * GENE_LENGTH

    # Population and time parameters
    POPULATION_SIZE = 500
    GENERATIONS = 100

    # Define mutation rates. Hypermutators have a significantly higher rate.
    # This represents the chance for a single base to mutate per generation.
    WT_MUTATION_RATE = 0.0005
    HYPERMUTATOR_MUTATION_RATE = 0.01  # 20x higher than wild-type

    def simulate_population(mutation_rate):
        """Helper function to run the simulation for one population."""
        total_mucoid_variants = 0
        unique_mucoid_genotypes = set()

        for _ in range(GENERATIONS):
            for _ in range(POPULATION_SIZE):
                # Each bacterium has a chance to mutate its mucA gene
                gene_list = list(WILD_TYPE_GENE)
                mutated = False
                for i in range(GENE_LENGTH):
                    if random.random() < mutation_rate:
                        gene_list[i] = 'B'  # 'B' represents a mutation
                        mutated = True
                
                if mutated:
                    total_mucoid_variants += 1
                    unique_mucoid_genotypes.add("".join(gene_list))
        
        return total_mucoid_variants, len(unique_mucoid_genotypes)

    print("--- Simulating Bacterial Adaptation ---")
    
    # Run simulations for both populations
    wt_frequency, wt_spectrum = simulate_population(WT_MUTATION_RATE)
    hm_frequency, hm_spectrum = simulate_population(HYPERMUTATOR_MUTATION_RATE)

    print("\nSimulation Results:")
    print("A 'mucoid variant' is a bacterium with at least one mutation in the mucA gene.")
    print("-" * 50)
    print("Wild-Type Population:")
    print(f"  - Total mucoid variants found (Frequency): {wt_frequency}")
    print(f"  - Unique mucoid genotypes (Spectrum): {wt_spectrum}")
    print("\nHypermutator Population:")
    print(f"  - Total mucoid variants found (Frequency): {hm_frequency}")
    print(f"  - Unique mucoid genotypes (Spectrum): {hm_spectrum}")
    print("-" * 50)
    
    print("\nConclusion from Simulation:")
    print("The final comparison shows:")
    print(f"Frequency: The hypermutator population produced a higher frequency of mucoid variants ({hm_frequency} > {wt_frequency}).")
    print(f"Spectrum: The hypermutator population produced a wider spectrum of mutations ({hm_spectrum} > {wt_spectrum}).")
    print("\nThis supports the choice that hypermutator strains lead to an increased frequency of mucoid variants with a wider spectrum of mucA mutations.")

# Execute the simulation
run_bacterial_evolution_simulation()
<<<B>>>
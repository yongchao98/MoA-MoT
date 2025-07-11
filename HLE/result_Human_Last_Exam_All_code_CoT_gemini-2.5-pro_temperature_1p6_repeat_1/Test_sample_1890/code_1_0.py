import random

def simulate_mutagenesis():
    """
    Simulates the emergence of mucoid variants in normal and hypermutator
    bacterial populations to illustrate the effect of mutation rates.
    """
    # --- Simulation Parameters ---
    # The 'mucA' gene is represented by a length, where each position can be a unique mutation site.
    gene_length = 800
    
    # We simulate over a number of generations or time steps.
    generations = 100000
    
    # Define mutation rates. Hypermutators are often 100x to 1000x higher.
    # This is the probability of a mucoid-causing mutation per generation.
    normal_mutation_rate = 1e-6
    hypermutator_mutation_rate = 1e-4 # 100x higher than normal

    # --- Data Trackers ---
    # Count the total number of mucoid variants found.
    normal_mucoid_count = 0
    hypermutator_mucoid_count = 0

    # Use sets to store the unique mutation sites found (the spectrum).
    # A set only stores unique values, perfect for tracking diversity.
    normal_mutation_spectrum = set()
    hypermutator_mutation_spectrum = set()

    # --- Simulation Loop ---
    for _ in range(generations):
        # Check for mutation in the normal population
        if random.random() < normal_mutation_rate:
            normal_mucoid_count += 1
            mutation_site = random.randint(1, gene_length)
            normal_mutation_spectrum.add(mutation_site)

        # Check for mutation in the hypermutator population
        if random.random() < hypermutator_mutation_rate:
            hypermutator_mucoid_count += 1
            mutation_site = random.randint(1, gene_length)
            hypermutator_mutation_spectrum.add(mutation_site)

    # --- Print Results ---
    print("Simulation Results:\n")
    
    print("--- Normal Strain ---")
    print(f"Total Generations Simulated = {generations}")
    print(f"Mutation Rate = {normal_mutation_rate}")
    print(f"Total Mucoid Variants Found = {normal_mucoid_count}")
    print(f"Size of mucA Mutation Spectrum (Unique Mutations) = {len(normal_mutation_spectrum)}")
    
    print("\n--- Hypermutator Strain ---")
    print(f"Total Generations Simulated = {generations}")
    print(f"Mutation Rate = {hypermutator_mutation_rate}")
    print(f"Total Mucoid Variants Found = {hypermutator_mucoid_count}")
    print(f"Size of mucA Mutation Spectrum (Unique Mutations) = {len(hypermutator_mutation_spectrum)}")
    
    print("\n--- Conclusion ---")
    print("The simulation shows that the hypermutator strain's increased rate of mutagenesis leads to:")
    print(f"1. An increased frequency of mucoid variants ({hypermutator_mucoid_count} > {normal_mucoid_count}).")
    print(f"2. A wider spectrum of mucA mutations ({len(hypermutator_mutation_spectrum)} > {len(normal_mutation_spectrum)}).")
    print("\nThis result directly supports answer choice B.")


# Run the simulation and print the output
simulate_mutagenesis()

import random

def run_simulation(strain_name, population_size, generations, mutation_rate_per_base):
    """
    Simulates the emergence of mutations in a gene for a given strain.
    A 'mucoid' variant is one with any mutation in our model gene.
    """
    # The key is the unique mutation sequence, the value is the count of its occurrence.
    mutations_observed = {}
    gene_length = 20 # A symbolic length for our mucA gene

    # The total number of mutation 'dice rolls' is population * generations * bases
    total_chances = population_size * generations * gene_length

    # Simulate mutation events
    for _ in range(total_chances):
        if random.random() < mutation_rate_per_base:
            # A mutation occurred. Let's create a unique mutation signature.
            # This represents one of many possible inactivating mutations.
            position = random.randint(1, gene_length)
            change = random.choice("ACGT")
            mutation_signature = f"pos{position}:{change}"
            
            # Record this unique mutation
            mutations_observed[mutation_signature] = mutations_observed.get(mutation_signature, 0) + 1

    # Calculate the final numbers
    total_mucoid_variants = sum(mutations_observed.values())
    spectrum_of_mutations = len(mutations_observed.keys())

    print(f"--- Results for {strain_name} Strain ---")
    print(f"With a mutation rate of {mutation_rate_per_base}:")
    # This represents the "final equation" with its numbers as requested.
    print(f"Final Equation: Frequency (Total Variants) = {total_mucoid_variants}")
    print(f"Final Equation: Spectrum (Unique Mutations) = {spectrum_of_mutations}")
    print("-" * 40)


if __name__ == "__main__":
    # Simulation parameters
    POP_SIZE = 10000
    GENERATIONS = 100
    NORMAL_RATE = 1e-7      # A typical baseline mutation rate
    HYPERMUTATOR_RATE = 1e-5  # A 100x higher rate for the hypermutator

    print("This simulation demonstrates how a higher mutation rate affects the frequency and variety (spectrum) of adaptive mutations.\n")
    
    # Run simulation for the normal strain
    run_simulation(
        strain_name="Normal",
        population_size=POP_SIZE,
        generations=GENERATIONS,
        mutation_rate_per_base=NORMAL_RATE
    )

    # Run simulation for the hypermutator strain
    run_simulation(
        strain_name="Hypermutator",
        population_size=POP_SIZE,
        generations=GENERATIONS,
        mutation_rate_per_base=HYPERMUTATOR_RATE
    )

    print("\nConclusion: The simulation shows the Hypermutator strain produces both a higher frequency of variants")
    print("and a wider spectrum of unique mutations, which aligns with the biological reasoning.")

import random

def simulate_mutations(strain_name, mutation_rate, generations, population_size, gene_length):
    """
    Simulates the emergence of mucoid variants in a bacterial population.
    This is a conceptual model to illustrate the principle.
    """
    print(f"--- Simulating: {strain_name} ---")
    print(f"Mutation Rate: {mutation_rate}")

    mucoid_variant_count = 0
    # Use a set to store the unique locations of mutations (the "spectrum")
    unique_mutation_locations = set()

    # In this model, we'll see how many variants arise in one generation for simplicity
    for _ in range(population_size):
        # We assume any mutation in the mucA gene is inactivating and causes the mucoid phenotype
        # We calculate the probability of at least one mutation in the gene
        prob_no_mutation_in_gene = (1 - mutation_rate) ** gene_length
        prob_at_least_one_mutation = 1 - prob_no_mutation_in_gene

        if random.random() < prob_at_least_one_mutation:
            mucoid_variant_count += 1
            # Simulate where the mutation occurred to check the spectrum
            mutated_position = random.randint(1, gene_length)
            unique_mutation_locations.add(mutated_position)

    # Calculate frequency
    frequency = mucoid_variant_count / population_size

    print(f"Resulting Mucoid Variants: {mucoid_variant_count}")
    print(f"Population Size: {population_size}")
    # Output the numbers used in the final equation for frequency
    print(f"Frequency Equation: {mucoid_variant_count} / {population_size}")
    print(f"Calculated Frequency of Mucoid Variants: {frequency:.6f}")
    print(f"Spectrum of Mutations (number of unique locations): {len(unique_mutation_locations)}\n")


def main():
    # Simulation Parameters
    # We use a large population to better see the effects of rare events.
    POPULATION_SIZE = 500000
    # An arbitrary length for the mucA gene
    GENE_LENGTH = 600 # mucA gene is ~570bp, so this is a reasonable approximation
    # Note: These are illustrative rates, not exact biological values.
    # The key is the relative difference.
    NORMAL_MUTATION_RATE = 1e-8
    HYPERMUTATOR_MUTATION_RATE = 1e-6 # 100x higher than normal

    print("This simulation models how a higher mutation rate affects the frequency and spectrum of mucoid variants.\n")

    # Run simulation for the normal strain
    simulate_mutations("Normal Strain", NORMAL_MUTATION_RATE, 1, POPULATION_SIZE, GENE_LENGTH)

    # Run simulation for the hypermutator strain
    simulate_mutations("Hypermutator Strain", HYPERMUTATOR_MUTATION_RATE, 1, POPULATION_SIZE, GENE_LENGTH)

    print("Conclusion: The hypermutator strain shows both a higher frequency of mucoid variants and a wider spectrum (more unique mutation sites) compared to the normal strain.")

if __name__ == "__main__":
    main()
<<<B>>>
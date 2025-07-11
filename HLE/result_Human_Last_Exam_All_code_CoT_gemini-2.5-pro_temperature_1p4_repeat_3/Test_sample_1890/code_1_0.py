import random

def run_bacterial_evolution_simulation(mutation_rate, generations, population_size):
    """
    Simulates the emergence of mucoid variants based on a given mutation rate.
    
    Returns:
        - The total number of mucoid variants generated.
        - The number of unique mutations (the "spectrum").
    """
    mucoid_variant_count = 0
    # Use a set to store unique mutations
    unique_mucA_mutations = set()

    for _ in range(generations):
        for _ in range(population_size):
            # Check if a random mutation occurs in this bacterium
            if random.random() < mutation_rate:
                # Assume the mutation causes a mucoid phenotype by affecting mucA
                mucoid_variant_count += 1
                
                # Generate a random identifier for the specific mutation site/type
                # to simulate a spectrum of different possible mutations.
                mutation_id = random.randint(1, 100000)
                unique_mucA_mutations.add(mutation_id)
                
    return mucoid_variant_count, len(unique_mucA_mutations)

# Simulation parameters
generations = 500
population_size = 10000

# Define mutation rates
# A hypermutator strain might have a 100x higher mutation rate
normal_rate = 0.00001
hypermutator_rate = 0.001

# Run simulation for the normal strain
normal_variants, normal_spectrum = run_bacterial_evolution_simulation(normal_rate, generations, population_size)

# Run simulation for the hypermutator strain
hypermutator_variants, hypermutator_spectrum = run_bacterial_evolution_simulation(hypermutator_rate, generations, population_size)

print("--- Conceptual Simulation Results ---")
print(f"Normal Strain:")
print(f"  - Frequency (Count): {normal_variants} mucoid variants")
print(f"  - Spectrum (Unique Mutations): {normal_spectrum} types")
print("\n")
print(f"Hypermutator Strain:")
print(f"  - Frequency (Count): {hypermutator_variants} mucoid variants")
print(f"  - Spectrum (Unique Mutations): {hypermutator_spectrum} types")
print("\n--- Conclusion ---")
print("The simulation shows that the hypermutator strain's increased mutation rate leads to an...")
print("-> Increased frequency of mucoid variants AND")
print("-> A wider spectrum of mucA mutations.")
print("This corresponds to answer choice B.")

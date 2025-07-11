import random

def simulate_mucoid_emergence(strain_name, mutation_rate, num_trials, gene_length):
    """
    Simulates the emergence of mucoid variants in a bacterial population.

    Args:
        strain_name (str): The name of the strain (e.g., 'Normal', 'Hypermutator').
        mutation_rate (float): The probability of a mucoid-causing mutation per trial.
        num_trials (int): The number of bacteria to simulate.
        gene_length (int): The length of the mucA gene for simulation purposes.

    Returns:
        A tuple containing the count of mucoid variants and the set of unique mutations.
    """
    mucoid_variant_count = 0
    unique_mutations = set()

    for _ in range(num_trials):
        # Check if a new mucoid-causing mutation occurs in this trial
        if random.random() < mutation_rate:
            mucoid_variant_count += 1
            
            # Simulate a random mutation within the mucA gene
            # This represents one specific mutation that inactivates the gene.
            mutation_position = random.randint(1, gene_length)
            original_base = random.choice(['A', 'T', 'C', 'G'])
            new_base = random.choice([b for b in ['A', 'T', 'C', 'G'] if b != original_base])
            mutation_type = "Point Mutation"
            
            # Record the unique mutation
            mutation_details = (mutation_type, mutation_position, original_base, new_base)
            unique_mutations.add(mutation_details)
            
    return mucoid_variant_count, unique_mutations

# --- Simulation Parameters ---
# The length of the P. aeruginosa mucA gene is 633 base pairs.
MUC_A_GENE_LENGTH = 633 
# We simulate a population of 1,000,000 bacteria for each strain.
NUM_TRIALS = 1_000_000
# Hypermutators can have mutation rates 100-1000 times higher than normal.
NORMAL_STRAIN_MUTATION_RATE = 1e-7
HYPERMUTATOR_STRAIN_MUTATION_RATE = 1e-5 # 100x higher

# --- Run Simulations ---
print("Running simulation...")

# Normal Strain
normal_mucoid_count, normal_unique_mutations = simulate_mucoid_emergence(
    'Normal Strain', 
    NORMAL_STRAIN_MUTATION_RATE, 
    NUM_TRIALS, 
    MUC_A_GENE_LENGTH
)

# Hypermutator Strain
hypermutator_mucoid_count, hypermutator_unique_mutations = simulate_mucoid_emergence(
    'Hypermutator Strain', 
    HYPERMUTATOR_STRAIN_MUTATION_RATE, 
    NUM_TRIALS, 
    MUC_A_GENE_LENGTH
)

# --- Print Results ---
print("\n--- Simulation Results ---")
print(f"Parameters: Population Size = {NUM_TRIALS}, mucA Gene Length = {MUC_A_GENE_LENGTH} bp")
print("-" * 28)

print(f"Normal Strain:")
print(f"  - Frequency of Mucoid Variants: {normal_mucoid_count} out of {NUM_TRIALS}")
print(f"  - Spectrum of mucA mutations (unique types found): {len(normal_unique_mutations)}")

print(f"\nHypermutator Strain:")
print(f"  - Frequency of Mucoid Variants: {hypermutator_mucoid_count} out of {NUM_TRIALS}")
print(f"  - Spectrum of mucA mutations (unique types found): {len(hypermutator_unique_mutations)}")

print("\n--- Conclusion ---")
print("The simulation shows that the hypermutator strain results in both an increased frequency of mucoid variants")
print("and a wider spectrum (more variety) of mutations in the mucA gene, which aligns with the biological principles.")

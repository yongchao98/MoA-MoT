import random

def simulate_mutations(mutation_rate, generations, gene_length):
    """
    Simulates mutations in a gene over a number of generations.
    Returns the total number of mutation events and a dictionary of mutation locations.
    """
    mutations_found = {}
    total_events = 0
    for _ in range(generations):
        # Check if a mutation occurs in this generation
        if random.random() < mutation_rate:
            total_events += 1
            # Determine where in the gene the mutation occurs
            position = random.randint(1, gene_length)
            mutations_found[position] = mutations_found.get(position, 0) + 1
    return total_events, mutations_found

# --- Parameters ---
# Hypermutator strains can have a 100- to 1000-fold increased mutation rate.
NORMAL_RATE = 0.0001
HYPERMUTATOR_RATE = 0.01  # 100x higher than normal
GENERATIONS = 50000
MUC_A_GENE_LENGTH = 573 # Approximate length of mucA gene in base pairs

print("Analyzing the effect of a hypermutator strain on the mucA gene...")
print(f"Simulation Parameters:\n- Normal Mutation Rate: {NORMAL_RATE}\n- Hypermutator Rate: {HYPERMUTATOR_RATE}\n- Generations: {GENERATIONS}\n")

# --- Simulation ---
normal_freq, normal_mutations = simulate_mutations(NORMAL_RATE, GENERATIONS, MUC_A_GENE_LENGTH)
hypermutator_freq, hypermutator_mutations = simulate_mutations(HYPERMUTATOR_RATE, GENERATIONS, MUC_A_GENE_LENGTH)

normal_spectrum = len(normal_mutations)
hypermutator_spectrum = len(hypermutator_mutations)

# --- Output Results ---
print("--- Results for Normal Strain ---")
# The frequency of mucoid variants is proportional to the number of mutation events.
print(f"Equation for Frequency: Total Events = {normal_freq}")
# The spectrum of mucA mutations is the number of unique sites mutated.
print(f"Equation for Spectrum: Unique Mutation Sites = {normal_spectrum}\n")


print("--- Results for Hypermutator Strain ---")
# The frequency is higher in the hypermutator strain.
print(f"Equation for Frequency: Total Events = {hypermutator_freq}")
# The spectrum is wider in the hypermutator strain.
print(f"Equation for Spectrum: Unique Mutation Sites = {hypermutator_spectrum}\n")

print("Conclusion:")
print("The hypermutator strain leads to a higher frequency of mucoid variants (more total events) and a wider spectrum of mucA mutations (more unique mutation sites).")

<<<B>>>
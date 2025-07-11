import random

# --- Parameters ---
# The probability that a single gene copy acquires a loss-of-function mutation in one event.
# Let's set it to a relatively high value for a clear demonstration.
p_mutation = 0.05

# The number of redundant copies in the multigene family.
n_copies_in_family = 4

# Number of simulation runs to get a stable average.
n_simulations = 10000

# --- Function to simulate the fate of a gene ---
def simulate_gene_loss(n_copies, p_mutation):
    """
    Simulates if a gene function is lost.
    A function is lost only if ALL copies are mutated.
    Returns True if function is lost, False otherwise.
    """
    for _ in range(n_copies):
        # If we find at least one copy that did NOT mutate, the function is preserved.
        if random.random() > p_mutation:
            return False
    # If we went through all copies and all of them mutated, the function is lost.
    return True

# --- Main Simulation ---
loss_count_single_gene = 0
loss_count_multigene_family = 0

for _ in range(n_simulations):
    # Scenario 1: Single-copy gene
    if simulate_gene_loss(1, p_mutation):
        loss_count_single_gene += 1

    # Scenario 2: Multigene family
    if simulate_gene_loss(n_copies_in_family, p_mutation):
        loss_count_multigene_family += 1

# --- Output Results ---
print("--- Simulation of Functional Gene Loss ---")
print(f"Parameters:")
print(f"  Probability of a single gene mutating (p): {p_mutation}")
print(f"  Number of copies in multigene family (n): {n_copies_in_family}")
print(f"  Number of simulations run: {n_simulations}\n")

# This section prints the numbers in the final theoretical equation.
print("--- Theoretical Probability (Equation: P(loss) = p^n) ---")
# Equation for the single-copy gene
p_loss_single_theoretical = p_mutation ** 1
print(f"Single-copy gene (n=1):")
print(f"  Equation: {p_mutation:.2f} ^ 1 = {p_loss_single_theoretical:.6f}")

# Equation for the multigene family
p_loss_family_theoretical = p_mutation ** n_copies_in_family
print(f"Multigene family (n={n_copies_in_family}):")
print(f"  Equation: {p_mutation:.2f} ^ {n_copies_in_family} = {p_loss_family_theoretical:.6f}\n")

print("--- Simulation Results ---")
# Calculate the observed probability (frequency) from the simulation
p_loss_single_simulated = loss_count_single_gene / n_simulations
p_loss_multigene_simulated = loss_count_multigene_family / n_simulations

print(f"Observed loss frequency for single-copy gene: {p_loss_single_simulated:.6f} ({loss_count_single_gene}/{n_simulations})")
print(f"Observed loss frequency for multigene family: {p_loss_multigene_simulated:.6f} ({loss_count_multigene_family}/{n_simulations})\n")

print("Conclusion: The multigene family provides significant protection against loss of function,")
print("compensating for the inability to remove deleterious mutations via recombination.")
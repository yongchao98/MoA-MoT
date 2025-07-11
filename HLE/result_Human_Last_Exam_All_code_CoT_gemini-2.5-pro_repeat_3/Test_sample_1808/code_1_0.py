import sys
# Redirect print to a string to capture output for the final answer format
# This is a helper function to format the output as requested.
# In a real scenario, you would just use regular print().
original_stdout = sys.stdout
from io import StringIO
string_io = StringIO()
sys.stdout = string_io

def calculate_fst(p1, p2):
    """Calculates Fst for a two-population system."""
    if p1 < 0 or p1 > 1 or p2 < 0 or p2 > 1:
        return float('nan')
    q1 = 1 - p1
    q2 = 1 - p2

    # Expected heterozygosity within each subpopulation
    hs1 = 2 * p1 * q1
    hs2 = 2 * p2 * q2

    # Average expected heterozygosity across subpopulations
    hs_avg = (hs1 + hs2) / 2

    # Expected heterozygosity in the total (metapopulation)
    p_total = (p1 + p2) / 2
    q_total = 1 - p_total
    ht = 2 * p_total * q_total

    # Fst calculation
    if ht == 0:
        return 0 # No variation in the total population
    fst = (ht - hs_avg) / ht
    return fst, ht, hs_avg

# --- Simulation Parameters ---
# Initial allele 'p' frequency in population 1
p1 = 1.0
# Initial allele 'p' frequency in population 2
p2 = 0.0
# Migration rate (proportion of a population replaced by migrants each generation)
m = 0.05
# Number of generations to simulate
generations = 50

# --- Simulation Start ---
print("--- Simulation of Gene Flow's Effect on Fst ---")
print(f"Initial conditions: Pop1 p={p1}, Pop2 p={p2}, Migration rate m={m}\n")

# Calculate initial Fst
fst, ht, hs = calculate_fst(p1, p2)
print(f"Generation 0:")
print(f"  Allele Frequencies: p1={p1:.4f}, p2={p2:.4f}")
print(f"  Fst Calculation: Ht={ht:.4f}, Hs={hs:.4f}")
print(f"  Fst = ({ht:.4f} - {hs:.4f}) / {ht:.4f} = {fst:.4f}\n")


# Run simulation for a few generations
for gen in range(1, generations + 1):
    # Allele frequencies change due to migration
    p1_new = (1 - m) * p1 + m * p2
    p2_new = (1 - m) * p2 + m * p1
    p1, p2 = p1_new, p2_new

    # Print status at intervals
    if gen in [1, 5, 10, 25, 50]:
        fst, ht, hs = calculate_fst(p1, p2)
        print(f"Generation {gen}:")
        print(f"  Allele Frequencies: p1={p1:.4f}, p2={p2:.4f}")
        # For the final generation, print the full equation as requested
        if gen == 50:
            print(f"  Final Fst Calculation: Ht={ht:.4f}, Hs={hs:.4f}")
            print(f"  Final Fst = ({ht:.4f} - {hs:.4f}) / {ht:.4f} = {fst:.4f}\n")
        else:
             print(f"  Fst = {fst:.4f}\n")

print("Conclusion: As shown, gene flow (m > 0) consistently reduces Fst over time.")
print("Therefore, high gene flow and high Fst are opposing forces.")

# Restore original stdout and print captured output
sys.stdout = original_stdout
captured_output = string_io.getvalue()
print(captured_output)
<<<A>>>
# Plan:
# 1. Define hypothetical fitness values (e.g., growth rates) for parental lines, 
#    F1 hybrids (from initial gene flow), and F2 hybrids (after meiosis).
# 2. Set the parental fitness as the reference (relative fitness = 1.0).
# 3. Calculate the relative fitness (w) for F1 and F2 hybrids.
# 4. Calculate the selection coefficient (s) for both generations, which quantifies the cost of gene flow.
# 5. Print the results clearly, showing the equation for each step.

# Step 1: Define hypothetical fitness values
parental_fitness = 1.0  # Fitness of the non-hybrid parent line (our reference)
f1_hybrid_fitness = 0.95 # Fitness of the first-generation hybrid
f2_hybrid_fitness = 0.70 # Fitness of the second-generation hybrid (after meiosis)

print("This script calculates the selection coefficient (s) to measure the cost of gene flow.\n")
print(f"Assumed Fitness Values (e.g., relative growth rate):")
print(f"  - Parental Lines: {parental_fitness}")
print(f"  - F1 Hybrids: {f1_hybrid_fitness}")
print(f"  - F2 Hybrids (post-meiosis): {f2_hybrid_fitness}\n")

# Step 2 & 3: Calculate relative fitness (w) for F1 and F2 hybrids
# Relative fitness w = (absolute fitness of genotype) / (absolute fitness of reference genotype)
w_f1 = f1_hybrid_fitness / parental_fitness
w_f2 = f2_hybrid_fitness / parental_fitness

print("--- Calculating Cost in F1 Generation ---")
print(f"Relative Fitness (w_f1) = {f1_hybrid_fitness} / {parental_fitness} = {w_f1:.2f}")

# Step 4 & 5: Calculate selection coefficient (s = 1 - w) and print
s_f1 = 1 - w_f1
print(f"Selection Coefficient (s_f1) = 1 - {w_f1:.2f} = {s_f1:.2f}")
print(f"This represents a {s_f1:.0%} fitness cost in the F1 generation.\n")

print("--- Calculating Cost in F2 Generation (after 'within mating' and meiosis) ---")
print(f"Relative Fitness (w_f2) = {f2_hybrid_fitness} / {parental_fitness} = {w_f2:.2f}")

s_f2 = 1 - w_f2
print(f"Selection Coefficient (s_f2) = 1 - {w_f2:.2f} = {s_f2:.2f}")
print(f"This represents a {s_f2:.0%} fitness cost in the F2 generation, revealing the effects of meiosis.")
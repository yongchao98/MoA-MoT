import math

# This script demonstrates how to calculate the cost of gene flow (selection coefficient, s)
# based on fitness measurements, as described in the best answer choice.

# --- Hypothetical Fitness Data ---
# Fitness can be measured in various ways, e.g., growth rate, biomass yield, etc.
# Here, we use a single value to represent the overall fitness of a population.
# A fitness of 1.0 is the reference (the non-hybridized parent).

# Fitness of the non-hybridized parental line (control group with no gene flow)
fitness_parent = 1.0

# Fitness of the first-generation (F1) hybrid.
# Sometimes F1 hybrids show 'hybrid vigor' and have fitness equal to or greater than the parent.
fitness_f1_hybrid = 1.02

# Fitness of the second-generation (F2) hybrid, after 'within mating' of F1s.
# This is where outbreeding depression (cost of gene flow) often appears due to meiosis.
fitness_f2_hybrid = 0.85

print("--- Calculating the Cost of Gene Flow (Selection Coefficient 's') ---\n")

# --- Calculation for F1 Hybrids ---
# The relative fitness (W) of the hybrid is its fitness compared to the parent.
relative_fitness_f1 = fitness_f1_hybrid / fitness_parent
# The selection coefficient (s) is calculated as 1 - W.
# A negative 's' indicates a fitness advantage (heterosis/hybrid vigor).
selection_coefficient_f1 = 1 - relative_fitness_f1

print("Step 1: Analyzing the F1 Hybrid Generation")
print(f"Parent Fitness: {fitness_parent}")
print(f"F1 Hybrid Fitness: {fitness_f1_hybrid}")
print(f"Calculating Relative Fitness (W) of F1 Hybrid: W = {fitness_f1_hybrid} / {fitness_parent} = {relative_fitness_f1:.2f}")
print(f"Calculating Selection Coefficient (s) for F1: s = 1 - {relative_fitness_f1:.2f} = {selection_coefficient_f1:.2f}")
if selection_coefficient_f1 < 0:
    print("Result: F1 generation shows a fitness advantage (hybrid vigor).\n")
else:
    print(f"Result: F1 generation shows a fitness cost of {selection_coefficient_f1:.2%}.\n")


# --- Calculation for F2 Hybrids (Post-Meiosis) ---
# We do the same calculation for the F2 generation to see the effects of meiotic recombination.
relative_fitness_f2 = fitness_f2_hybrid / fitness_parent
# A positive 's' indicates a fitness disadvantage (cost/outbreeding depression).
selection_coefficient_f2 = 1 - relative_fitness_f2

print("Step 2: Analyzing the F2 Hybrid Generation (after meiosis)")
print(f"Parent Fitness: {fitness_parent}")
print(f"F2 Hybrid Fitness: {fitness_f2_hybrid}")
print(f"Calculating Relative Fitness (W) of F2 Hybrid: W = {fitness_f2_hybrid} / {fitness_parent} = {relative_fitness_f2:.2f}")
print(f"Calculating Selection Coefficient (s) for F2: s = 1 - {relative_fitness_f2:.2f} = {selection_coefficient_f2:.2f}")
if selection_coefficient_f2 > 0:
    print(f"Result: F2 generation shows a fitness cost (outbreeding depression) of {selection_coefficient_f2:.2%}.\n")
else:
    print("Result: F2 generation shows no fitness cost.\n")

print("Conclusion: The cost of gene flow is most clearly observed in the F2 generation after meiotic recombination.")

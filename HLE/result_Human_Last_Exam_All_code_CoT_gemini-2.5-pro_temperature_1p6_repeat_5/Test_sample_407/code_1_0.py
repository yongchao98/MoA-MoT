# Plan:
# 1. Define hypothetical fitness values for a parental line (no gene flow)
#    and a hybrid F2 line (after gene flow and subsequent meiosis).
#    Fitness can be measured by growth rate, for example.
# 2. Define the formula for the selection coefficient (s), which is a measure
#    of the cost. The formula is s = 1 - (W_hybrid / W_parental), where W is fitness.
# 3. Calculate the selection coefficient using the hypothetical values.
# 4. Print the entire calculation step-by-step to show how the cost is quantified.

# Step 1: Define hypothetical fitness values.
# Let's assume fitness is measured as the exponential growth rate of the yeast colony.
# W_parental represents the average fitness of the non-hybrid control lines.
W_parental = 0.55  # e.g., in units of doublings per hour

# W_hybrid_F2 represents the average fitness of the F2 generation of the hybrids.
# We expect this to be lower if there is a cost due to hybrid breakdown after meiosis.
W_hybrid_F2 = 0.48 # e.g., in units of doublings per hour

# Step 2: The selection coefficient 's' quantifies the fitness cost.
# It represents the selection pressure against the hybrids.
# Formula: s = 1 - (relative fitness of hybrids)

# Step 3 & 4: Calculate and print the result.
relative_fitness = W_hybrid_F2 / W_parental
selection_coefficient = 1 - relative_fitness

print("Measuring the Cost of Gene Flow (Selection Coefficient)")
print("-----------------------------------------------------")
print(f"Fitness of parental (no gene flow) line (W_parental): {W_parental}")
print(f"Fitness of hybrid F2 line (W_hybrid_F2): {W_hybrid_F2}")
print("\nThe selection coefficient (s) is calculated as: s = 1 - (W_hybrid_F2 / W_parental)")
print("\nCalculation:")
print(f"s = 1 - ({W_hybrid_F2} / {W_parental})")
print(f"s = 1 - {relative_fitness:.4f}")
print(f"s = {selection_coefficient:.4f}")

print("\nA positive selection coefficient indicates a fitness cost for the hybrids.")
print(f"In this case, the cost of gene flow is {selection_coefficient:.2%}.")
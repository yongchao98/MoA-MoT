import math

# Step 1: Define the molecular properties based on the analysis of cyclopentanecarboxylic acid.

# Total number of non-hydrogen atoms (6 Carbons + 2 Oxygens)
N = 8

# Sizes of the symmetry equivalence classes (ni).
# The order does not matter for the final sum.
# {C_ring-alpha}, {C_ring-beta, C_ring-beta'}, {C_ring-gamma, C_ring-gamma'}, {C_carboxyl}, {O_carbonyl}, {O_hydroxyl}
ni_list = [1, 2, 2, 1, 1, 1]

# Step 2: Calculate the two main terms of the Böttcher formula.
# First term: N * log2(N)
term1 = N * math.log2(N)

# Second term: The sum of ni*log2(ni) for each class
sum_term_parts = []
for ni in ni_list:
    # math.log2(1) is 0, so 1*log2(1) is 0.
    if ni > 1:
        sum_term_parts.append(ni * math.log2(ni))
    else:
        sum_term_parts.append(0.0)
term2 = sum(sum_term_parts)

# Step 3: Calculate the final complexity index.
complexity = term1 - term2

# Step 4: Print the detailed calculation as requested.
print("Calculation of Böttcher Molecular Complexity for Cyclopentanecarboxylic Acid:\n")
print(f"Formula: C = N*log2(N) - sum(ni*log2(ni))")
print(f"N (non-hydrogen atoms) = {N}")
print(f"ni (sizes of symmetry classes) = {ni_list}\n")

# Build and print the equation with substituted numbers
equation_part1 = f"{N}*log2({N})"

# We re-sort the list for a prettier equation printout, though it's not necessary for the calculation.
sorted_ni_list = sorted(ni_list, reverse=True)
equation_part2_terms = [f"{ni}*log2({ni})" for ni in sorted_ni_list]
equation_part2 = " + ".join(equation_part2_terms)

print("Final Equation:")
# The equation string is constructed to show each number.
print(f"C = ({equation_part1}) - ({equation_part2})")

# Print the evaluated terms
sorted_sum_term_parts = sorted([p for p in sum_term_parts if p > 0], reverse=True)
# Add the zero terms back for clarity
while len(sorted_sum_term_parts) < len(ni_list):
    sorted_sum_term_parts.append(0.0)
    
print(f"C = ({term1}) - ({' + '.join([str(p) for p in sorted_sum_term_parts])})")
print(f"C = {term1} - {term2}")
print(f"C = {complexity}")

print(f"\nThe Böttcher Molecular Complexity is {complexity}.")
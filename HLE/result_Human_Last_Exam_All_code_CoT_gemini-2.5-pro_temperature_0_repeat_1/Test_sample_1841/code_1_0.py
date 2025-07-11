import sys

# Plan:
# 1. Define the sample data: a list where each number is the count of individuals for a species.
#    To get D=0, we need a sample with only one species.
# 2. Calculate N, the total number of individuals.
# 3. Calculate the numerator of the Simpson's Index fraction: sum of n(n-1) for each species.
# 4. Calculate the denominator: N(N-1).
# 5. Calculate the final Simpson's Diversity Index, D.
# 6. Print the explanation and the final equation with the calculated numbers.

# Let's assume the student's sample consisted of 45 individuals of a single bat species.
species_counts = [45]
species_names = ["Species A"]

# Calculate N, the total number of individuals
N = sum(species_counts)

# The numerator of the fraction is the sum of n(n-1) for each species
# Since there is only one species, n = N
sum_n_n_minus_1 = 0
for n in species_counts:
    # This check prevents calculation with n<2, which would be invalid.
    if n > 1:
        sum_n_n_minus_1 += n * (n - 1)

# The denominator of the fraction is N(N-1)
# We must handle the case where N < 2, as N(N-1) would be 0 or negative.
if N < 2:
    print("Cannot calculate Simpson's Index: The total number of individuals (N) must be 2 or more.")
    sys.exit()

N_N_minus_1 = N * (N - 1)

# Calculate Simpson's Diversity Index (D)
# D = 1 - (numerator / denominator)
simpson_diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)

print("--- Analysis of Simpson's Diversity Index ---")
print(f"Sample Data: The student found {N} bats, all belonging to one species.")
print("\nFormula: D = 1 - [ Σn(n-1) / N(N-1) ]")
print("\nCalculations:")
print(f"n = The count for the single species found = {species_counts[0]}")
print(f"N = Total individuals = {N}")
print(f"Σn(n-1) = {species_counts[0]} * ({species_counts[0]} - 1) = {sum_n_n_minus_1}")
print(f"N(N-1) = {N} * ({N} - 1) = {N_N_minus_1}")

print("\nFinal Equation with calculated values:")
# Using an f-string to format the output clearly
print(f"D = 1 - ( {sum_n_n_minus_1} / {N_N_minus_1} ) = {simpson_diversity_index}")

print("\nConclusion:")
print("The result D=0 is mathematically valid, as shown by the calculation.")
print("However, it is not ecologically valid because it contradicts the known fact that multiple bat species exist on the island. The sample was not representative of the ecosystem.")

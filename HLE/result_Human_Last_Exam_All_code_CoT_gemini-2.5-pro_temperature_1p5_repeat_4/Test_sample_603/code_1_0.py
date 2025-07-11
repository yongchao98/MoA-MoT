import numpy as np

# This script illustrates the concept of heterogeneity lost in bulk measurements.
# We simulate a small population of DNA molecules that, due to minor
# conformational differences or local environments, have slightly different
# melting temperatures (Tm). This represents the "heterogeneity".

# 1. Define the melting temperatures (in Celsius) for a set of individual molecules.
# This represents the true, heterogeneous population that a single-molecule
# experiment could potentially distinguish.
single_molecule_tms = [89.8, 90.1, 90.3, 90.5, 91.2]

# 2. In a bulk experiment (like calorimetry), we don't see the individual values.
# Instead, we measure the average property of the entire population.
# We calculate the mean Tm to represent the result of a bulk measurement.
bulk_tm = np.mean(single_molecule_tms)

# 3. Print the results to show the difference.
print("Illustrating how bulk measurements average out heterogeneity:")
print("-" * 55)

# The detailed information that is lost
print("Hypothetical single-molecule melting temperatures (hidden heterogeneity):")
# Output each number from the initial list
for tm in single_molecule_tms:
    print(f"{tm}°C")

print("\nEquation for the 'Bulk' measurement (Average):")
# Create the equation string and print it
equation_str = " + ".join(map(str, single_molecule_tms))
print(f"({equation_str}) / {len(single_molecule_tms)} = {bulk_tm:.2f}°C")

# The final, averaged result seen in the experiment
print("\nResult from a 'bulk' experiment (the average value):")
print(f"Tm = {bulk_tm:.2f}°C")

print("\nConclusion: The bulk experiment reports a single value, obscuring the underlying")
print("distribution of individual molecular behaviors.")

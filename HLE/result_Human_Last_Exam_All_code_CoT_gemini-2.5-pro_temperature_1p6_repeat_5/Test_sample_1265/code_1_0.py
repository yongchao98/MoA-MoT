import sys

# Define parameters for the simulation
# Let's assume a very low frequency allele for drug resistance
true_allele_frequency = 0.005  # 0.5% frequency
# A typical error rate for high-fidelity sequencing (e.g., Illumina)
sequencing_error_rate = 0.001  # 0.1% error rate

# The two sequencing depths from the problem
depth_1 = 40
depth_2 = 80

print(f"Goal: Detect a rare allele with a true frequency of {true_allele_frequency*100}%.")
print(f"Assumption: The sequencing technology has an error rate of {sequencing_error_rate*100}%.\n")
print("-" * 50)

# --- Calculation for 40X depth ---
print(f"Case 1: Sequencing at {depth_1}X depth\n")

# Expected number of reads showing the TRUE rare allele
expected_true_reads_1 = depth_1 * true_allele_frequency

# Expected number of reads showing the allele due to SEQUENCING ERROR
# This happens on reads that were supposed to be the reference/wild-type
expected_error_reads_1 = depth_1 * (1 - true_allele_frequency) * sequencing_error_rate

# The total number of times we expect to see the allele
total_observed_1 = expected_true_reads_1 + expected_error_reads_1

print("The equation for total observed alleles is:")
print("Total = (Depth * True_Frequency) + (Depth * (1 - True_Frequency) * Error_Rate)\n")

print(f"Expected number of reads with the true allele:")
print(f"Equation: {depth_1} * {true_allele_frequency} = {expected_true_reads_1:.2f}")

print(f"\nExpected number of reads showing the allele due to error:")
print(f"Equation: {depth_1} * (1 - {true_allele_frequency}) * {sequencing_error_rate} = {expected_error_reads_1:.4f}")

print(f"\nAt {depth_1}X depth, we expect to see the allele {total_observed_1:.2f} times.")
print(f"The signal (true reads) to noise (error reads) ratio is {expected_true_reads_1/expected_error_reads_1:.2f}:1")
print("This count is very low, making it hard to distinguish from random noise.")

print("-" * 50)

# --- Calculation for 80X depth ---
print(f"Case 2: Sequencing at {depth_2}X depth\n")

# Expected number of reads showing the TRUE rare allele
expected_true_reads_2 = depth_2 * true_allele_frequency

# Expected number of reads showing the allele due to SEQUENCING ERROR
expected_error_reads_2 = depth_2 * (1 - true_allele_frequency) * sequencing_error_rate

# The total number of times we expect to see the allele
total_observed_2 = expected_true_reads_2 + expected_error_reads_2

print(f"Expected number of reads with the true allele:")
print(f"Equation: {depth_2} * {true_allele_frequency} = {expected_true_reads_2:.2f}")

print(f"\nExpected number of reads showing the allele due to error:")
print(f"Equation: {depth_2} * (1 - {true_allele_frequency}) * {sequencing_error_rate} = {expected_error_reads_2:.4f}")

print(f"\nAt {depth_2}X depth, we expect to see the allele {total_observed_2:.2f} times.")
print(f"The signal (true reads) to noise (error reads) ratio is {expected_true_reads_2/expected_error_reads_2:.2f}:1")
print("By doubling the depth, we have increased the expected count, making the signal stronger and more statistically reliable.")

print("-" * 50)
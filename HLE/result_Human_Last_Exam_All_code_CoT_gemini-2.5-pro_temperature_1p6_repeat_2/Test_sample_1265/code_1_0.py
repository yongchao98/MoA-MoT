import numpy as np

# --- Parameters ---
true_rare_allele_frequency = 0.01  # A rare allele at 1% frequency
sequencing_error_rate = 0.002      # A realistic 0.2% sequencing error rate
depths_to_simulate = [40, 80]
rare_allele = 'T'
reference_allele = 'A'
# Other bases that an error can convert to
other_error_bases = ['C', 'G']

print(f"Simulating the detection of a rare allele ('{rare_allele}') with a true frequency of {true_rare_allele_frequency*100}%.")
print(f"Assuming a sequencing error rate of {sequencing_error_rate*100}%.")
print("-" * 60)

# --- Simulation ---
for depth in depths_to_simulate:
    print(f"\nSimulating at {depth}X sequencing depth:")

    # Generate the 'true' alleles we would get at this depth based on real frequency
    true_alleles_sampled = np.random.choice(
        [rare_allele, reference_allele],
        size=depth,
        p=[true_rare_allele_frequency, 1 - true_rare_allele_frequency]
    )

    # Now, simulate the sequencing process, which can introduce errors
    observed_reads = []
    for allele in true_alleles_sampled:
        # Check if a random sequencing error occurs for this read
        if np.random.random() < sequencing_error_rate:
            # An error occurred. The base is read incorrectly.
            if allele == reference_allele:
                # The reference base was misread as the rare base or another base
                observed_reads.append(np.random.choice([rare_allele] + other_error_bases))
            else: # The rare allele was misread
                observed_reads.append(np.random.choice([reference_allele] + other_error_bases))
        else:
            # No error occurred, the read is correct
            observed_reads.append(allele)

    # Count the final observed alleles in our sequencing data
    observed_rare_count = observed_reads.count(rare_allele)
    observed_ref_count = observed_reads.count(reference_allele)
    
    # Calculate observed frequency
    if depth > 0:
        observed_frequency = observed_rare_count / depth
    else:
        observed_frequency = 0

    print(f"Total reads (Depth): {depth}")
    print(f"Observed count of rare allele '{rare_allele}': {observed_rare_count}")
    print(f"Observed count of reference allele '{reference_allele}': {observed_ref_count}")
    
    # Per the instruction, showing the numbers in the final equation
    print(f"\nFinal Equation: Observed Frequency = Observed Rare Count / Total Depth")
    print(f"Calculation: {observed_frequency:.4f} = {observed_rare_count} / {depth}")

    # A simple conclusion for this simulation run
    # At very low counts, it's hard to distinguish from noise.
    if observed_rare_count > 1:
        print("Result: The rare allele was detected with some confidence.")
    else:
        print("Result: The rare allele was likely missed or is indistinguishable from sequencing noise.")
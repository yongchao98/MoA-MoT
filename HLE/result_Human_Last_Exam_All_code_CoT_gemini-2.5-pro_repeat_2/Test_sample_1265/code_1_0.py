def calculate_expected_reads(depth, frequency):
    """
    Calculates the expected number of reads that will show a rare allele
    given a specific sequencing depth and allele frequency.

    Args:
        depth (int): The average sequencing depth.
        frequency (float): The frequency of the rare allele in the population.
    """
    expected_reads = depth * frequency
    print(f"For a sequencing depth of {depth}X and an allele frequency of {frequency*100}%:")
    # The final code must output each number in the final equation
    print(f"The calculation is: {depth} * {frequency} = {expected_reads:.1f} expected reads.")
    
    if expected_reads < 1:
        print("Result: This is less than one read on average. The allele may be missed entirely or be indistinguishable from a single random sequencing error.\n")
    else:
        print(f"Result: With an expectation of {expected_reads:.1f} reads, the allele is more likely to be detected and confirmed as a real variant, not just noise.\n")


# --- Main ---
print("This script demonstrates why increasing sequencing depth is critical for detecting low-frequency alleles, as proposed in the correct strategy.\n")

# Define the frequency of our hypothetical rare, drug-resistant allele.
# Let's say it exists in 2% of the microbial cells.
allele_frequency = 0.02

# Case 1: Standard sequencing depth of 40X
depth_standard = 40
calculate_expected_reads(depth_standard, allele_frequency)

# Case 2: Increased sequencing depth of 80X
depth_increased = 80
calculate_expected_reads(depth_increased, allele_frequency)

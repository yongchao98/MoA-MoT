def calculate_expected_reads(allele_frequency, sequencing_depth):
    """
    Calculates the expected number of reads showing a rare allele.
    This helps demonstrate why higher sequencing depth is critical for detecting
    low-frequency variants.
    """
    print(f"Analyzing detection at {sequencing_depth}X depth...")
    print(f"Given a very low allele frequency of: {allele_frequency}")

    # The expected number of reads is the depth multiplied by the frequency.
    expected_reads = sequencing_depth * allele_frequency

    # Print out the equation and the result
    print("Equation: Expected Reads = Sequencing Depth * Allele Frequency")
    print(f"Calculation: {expected_reads:.2f} = {sequencing_depth} * {allele_frequency}")
    if expected_reads < 1:
        print("Result: On average, you would expect to see this allele on less than one read.")
        print("This makes confident detection extremely difficult or impossible.\n")
    elif expected_reads < 5:
        print(f"Result: On average, you would expect to see this allele on about {expected_reads:.2f} reads.")
        print("This is detectable, but requires high-quality data to distinguish from errors.\n")
    else:
        print(f"Result: On average, you would expect to see this allele on about {expected_reads:.2f} reads, which is more robust for detection.\n")

# Define the parameters for a very low frequency allele
low_allele_frequency = 0.005 # 0.5% frequency

# Case 1: Standard sequencing depth (e.g., 40X)
depth1 = 40
calculate_expected_reads(low_allele_frequency, depth1)

# Case 2: Increased sequencing depth (e.g., 80X), as suggested in the options
depth2 = 80
calculate_expected_reads(low_allele_frequency, depth2)

print("As the simulation shows, doubling the depth from 40X to 80X doubles the expected number of reads for the rare allele. This directly increases the chance of successful identification, highlighting a key principle in the correct answer.")
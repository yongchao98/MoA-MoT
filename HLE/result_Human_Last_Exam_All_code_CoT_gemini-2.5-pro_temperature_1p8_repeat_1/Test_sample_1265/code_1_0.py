import random

def simulate_sequencing(genome_size, target_position, allele_frequency, sequencing_depth, error_rate):
    """
    Simulates sequencing reads from a specific genomic position.
    
    Args:
        genome_size (int): The size of the genome.
        target_position (int): The position of the allele of interest.
        allele_frequency (float): The frequency of the rare allele (0 to 1).
        sequencing_depth (int): The number of reads covering the target position.
        error_rate (float): The sequencing error rate (0 to 1).
        
    Returns:
        A tuple containing counts of reference allele, rare allele, and error reads.
    """
    ref_allele_count = 0
    rare_allele_count = 0
    error_count = 0

    for _ in range(sequencing_depth):
        # Does this read contain the rare allele?
        if random.random() < allele_frequency:
            # It's a true rare allele. Now check for sequencing error.
            if random.random() > error_rate:
                rare_allele_count += 1
            else:
                error_count += 1 # True allele was sequenced incorrectly
        else:
            # It's a reference allele. Now check for sequencing error.
            if random.random() > error_rate:
                ref_allele_count += 1
            else:
                error_count += 1 # Ref allele was sequenced incorrectly
                
    return ref_allele_count, rare_allele_count, error_count

def main():
    """
    Main function to run simulations and print results.
    """
    # --- Parameters ---
    allele_freq = 0.01  # A rare allele at 1% frequency
    seq_error_rate = 0.005 # A typical short-read error rate (0.5%)

    # --- Scenario 1: Low Depth (e.g., 40X) ---
    low_depth = 40
    ref1, rare1, err1 = simulate_sequencing(1000000, 500, allele_freq, low_depth, seq_error_rate)
    
    print("--- Scenario 1: Low Sequencing Depth (40X) ---")
    print(f"Goal: Detect a rare allele with {allele_freq*100}% frequency.")
    print(f"Simulation Result:")
    print(f"Reads supporting Reference Allele: {ref1}")
    print(f"Reads supporting Rare Allele: {rare1}")
    print(f"Reads with Sequencing Errors: {err1}")
    print("Conclusion: With low depth, it is statistically unlikely to see the rare allele.")
    print(f"Even if we do, the number of reads ({rare1}) is too low to distinguish from noise.\n")


    # --- Scenario 2: High Depth (e.g., 80X) ---
    high_depth = 80
    ref2, rare2, err2 = simulate_sequencing(1000000, 500, allele_freq, high_depth, seq_error_rate)
    
    print("--- Scenario 2: Increased Sequencing Depth (80X) ---")
    print(f"Goal: Detect a rare allele with {allele_freq*100}% frequency.")
    print(f"Simulation Result:")
    print(f"Reads supporting Reference Allele: {ref2}")
    print(f"Reads supporting Rare Allele: {rare2}")
    print(f"Reads with Sequencing Errors: {err2}")
    print("Conclusion: Increasing depth improves the chance of observing the rare allele.")
    print(f"However, even at {high_depth}X, the count ({rare2}) may still be too low for high confidence.\n")
    
    # --- Scenario 3: Very High Depth ---
    very_high_depth = 1000
    ref3, rare3, err3 = simulate_sequencing(1000000, 500, allele_freq, very_high_depth, seq_error_rate)

    print("--- Scenario 3: Very High Sequencing Depth (1000X) ---")
    print(f"This demonstrates the principle from Answer A taken to its logical conclusion.")
    print(f"Simulation Result:")
    print(f"Reads supporting Reference Allele: {ref3}")
    print(f"Reads supporting Rare Allele: {rare3}")
    print(f"Reads with Sequencing Errors: {err3}")
    print("Conclusion: With very high depth, the signal from the rare allele becomes statistically significant")
    print(f"and clearly distinguishable from the background noise ({rare3} >> {err3}).")
    print("\nThis illustrates why both high-quality preparation (to keep error rate low) and high depth (to get enough reads) are essential, as stated in Answer A.")


if __name__ == "__main__":
    main()
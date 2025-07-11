import random

def simulate_sequencing(depth, true_frequency, error_rate):
    """
    Simulates sequencing a single nucleotide position to a specified depth.

    Args:
        depth (int): The sequencing depth (e.g., 40, 80, 500).
        true_frequency (float): The true frequency of the minor allele (e.g., 0.02 for 2%).
        error_rate (float): The probability of a sequencing error for a single base.
    """
    major_allele = 'A'
    minor_allele = 'T' # The allele of interest (e.g., drug resistance)
    
    alleles_observed = {major_allele: 0, minor_allele: 0, 'Error': 0}
    
    print(f"\n--- Simulating for {depth}X Depth ---")
    print(f"Goal: Detect minor allele '{minor_allele}' with true frequency of {true_frequency*100}%.")
    
    for _ in range(depth):
        # Step 1: Determine the true allele for this read based on population frequency
        if random.random() < true_frequency:
            true_base = minor_allele
        else:
            true_base = major_allele
            
        # Step 2: Simulate potential sequencing error
        if random.random() < error_rate:
            # An error occurred. Let's assume it flips the base.
            read_base = minor_allele if true_base == major_allele else major_allele
            alleles_observed['Error'] += 1
        else:
            # No error occurred. The read base is the true base.
            read_base = true_base

        # Step 3: Count the observed allele in our sequencing data
        if read_base == major_allele:
            alleles_observed[major_allele] += 1
        else:
            alleles_observed[minor_allele] += 1
            
    # Step 4: Report the results
    minor_allele_count = alleles_observed[minor_allele]
    major_allele_count = alleles_observed[major_allele]
    total_reads = minor_allele_count + major_allele_count
    
    # Check if we detected anything and calculate observed frequency
    if total_reads > 0 and minor_allele_count > 0:
        observed_freq = (minor_allele_count / total_reads) * 100
        print(f"Result: Detected! Observed frequency: {observed_freq:.2f}%")
    else:
        print("Result: Not detected!")

    # Print the "equation" of observed alleles
    print(f"Final Read Counts = {minor_allele_count} '{minor_allele}' reads + {major_allele_count} '{major_allele}' reads")


if __name__ == '__main__':
    # --- Parameters ---
    # A very low frequency allele, e.g., 1.5% of the population has it
    resistance_allele_frequency = 0.015 
    
    # Standard high-fidelity (e.g., Illumina) sequencing error rate
    sequencing_error_rate = 0.005 # 0.5% error rate
    
    # Depths to compare
    depths_to_test = [40, 80, 500]

    for depth in depths_to_test:
        simulate_sequencing(depth, resistance_allele_frequency, sequencing_error_rate)

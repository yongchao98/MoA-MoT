import random

def find_rare_allele_from_reads():
    """
    This script simulates the final analysis step from Option A:
    "compare these reads to the alleles of interest" after increasing sequencing depth.
    """
    # --- Step 1: Define Parameters ---
    # We assume a successful library prep and increased sequencing depth.
    total_sequencing_depth = 80  # Increased from 40X to 80X
    actual_resistance_frequency = 0.05  # A hypothetical low frequency of 5%

    reference_allele = 'C'
    # The allele we are looking for, which may confer drug resistance
    allele_of_interest = 'T'

    print("Plan: Identify a low-frequency resistance allele in a microbial population.")
    print(f"Simulating a sequencing depth of {total_sequencing_depth}X for the target gene.\n")

    # --- Step 2: Simulate the Read Data ---
    # In a real experiment, this data would come from a BAM/SAM file after alignment.
    # We calculate the expected number of reads showing the resistance allele.
    num_resistance_reads = int(total_sequencing_depth * actual_resistance_frequency)
    num_reference_reads = total_sequencing_depth - num_resistance_reads

    # Create a simulated "pileup" of bases from all 80 reads at this one position.
    simulated_reads = [reference_allele] * num_reference_reads + [allele_of_interest] * num_resistance_reads
    random.shuffle(simulated_reads)
    
    # --- Step 3: Analyze the Reads to Identify the Allele ---
    # This is the core "comparison" step. We count each allele's occurrences.
    print("Analyzing simulated reads to find the allele of interest...")
    count_of_interest = simulated_reads.count(allele_of_interest)
    count_of_reference = simulated_reads.count(reference_allele)
    total_reads_counted = len(simulated_reads)

    print(f"Found {count_of_reference} reads matching the reference allele ('{reference_allele}').")
    print(f"Found {count_of_interest} reads matching the allele of interest ('{allele_of_interest}').\n")

    # --- Step 4: Calculate Frequency and Report ---
    # With high depth, we can be more confident that this isn't just a sequencing error.
    if total_reads_counted > 0:
        observed_frequency = count_of_interest / total_reads_counted
        print("Final Calculation of the Allele's Frequency:")
        # The final equation with numbers is shown here, as requested.
        print(f"Frequency = {count_of_interest} (allele of interest) / {total_reads_counted} (total reads) = {observed_frequency:.4f}")
        
        if observed_frequency > 0.01: # A typical detection threshold
            print("\nResult: Very low frequency allele successfully identified.")
        else:
            print("\nResult: Allele of interest not detected or below detection threshold.")
    else:
        print("Result: No reads were found at the target location.")

# Run the simulation
find_rare_allele_from_reads()
import random

def run_simulation():
    """
    Demonstrates how increasing sequencing depth improves the detection of
    low-frequency alleles.
    """
    # --- Parameters ---
    genome_size = 20000    # The size of the microbial genome in base pairs.
    read_length = 150      # The length of a single sequencing read.
    allele_pos = 10000     # The specific position of the drug resistance allele.
    ref_base = 'A'         # The reference (wild-type) base at the position.
    alt_base = 'T'         # The alternate (resistance) base at the position.
    allele_freq = 0.015    # The true frequency of the resistance allele (1.5%).

    print("This script simulates sequencing to detect a rare drug-resistance allele.\n")
    print(f"Goal: Detect a '{alt_base}' allele at position {allele_pos} with a true frequency of {allele_freq:.1%}.\n")

    # --- Run Simulation for both depths ---
    simulate_at_depth(40, genome_size, read_length, allele_pos, ref_base, alt_base, allele_freq)
    simulate_at_depth(80, genome_size, read_length, allele_pos, ref_base, alt_base, allele_freq)

def simulate_at_depth(depth, genome_size, read_length, allele_pos, ref_base, alt_base, allele_freq):
    """
    A helper function to run the simulation for a specific depth.
    """
    print(f"--- Simulating Sequencing at {depth}X Depth ---")

    # This equation calculates the total number of reads needed.
    # num_reads = (genome_size * depth) / read_length
    num_reads = (genome_size * depth) // read_length
    print(f"Total reads to simulate: ({genome_size} * {depth}) / {read_length} = {num_reads}")

    ref_reads_count = 0
    alt_reads_count = 0

    # Simulate each read
    for _ in range(num_reads):
        start_pos = random.randint(0, genome_size - read_length)
        end_pos = start_pos + read_length

        # Check if this read covers the allele's position
        if start_pos <= allele_pos < end_pos:
            # The read covers the position. Decide if it carries the rare allele.
            if random.random() < allele_freq:
                alt_reads_count += 1
            else:
                ref_reads_count += 1

    total_coverage = ref_reads_count + alt_reads_count
    
    # Final equation for observed frequency
    # observed_freq = alt_reads_count / total_coverage
    if total_coverage > 0:
        observed_freq = alt_reads_count / total_coverage
    else:
        observed_freq = 0

    print("\nResults:")
    print(f"Reads with Reference Base ('{ref_base}'): {ref_reads_count}")
    print(f"Reads with Alternate Base ('{alt_base}'): {alt_reads_count}")
    print(f"Total Coverage at Position {allele_pos}: {total_coverage} reads")
    
    print(f"Observed Frequency Equation: {alt_reads_count} / {total_coverage} = {observed_freq:.3%}")

    # A typical variant caller requires at least 3-5 supporting reads to make a confident call.
    if alt_reads_count < 3:
        print("Conclusion: Detection is UNLIKELY. The evidence is too weak and could be mistaken for error.")
    else:
        print("Conclusion: Detection is LIKELY. There is sufficient evidence for the rare allele.")
    print("-" * 45 + "\n")

if __name__ == '__main__':
    run_simulation()
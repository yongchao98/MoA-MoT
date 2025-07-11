import random

def solve_recombination_frequency():
    """
    This function simulates bacterial conjugation to determine which genetic marker
    would have the highest recombination frequency.

    Plan:
    1. Define the genes and their order of transfer based on the problem (thr-azi-gal).
       We assign arbitrary time units for entry into the recipient cell.
    2. Simulate a large number of interrupted mating events. In each simulation,
       the duration of DNA transfer is random, as the connection between cells can
       break at any time.
    3. For each simulated mating, we count which genes were successfully transferred
       based on the random mating duration. A gene is transferred only if the mating
       duration is longer than its required entry time.
    4. Calculate the frequency of transfer for each gene. This frequency is directly
       proportional to the frequency of recombinants observed in an experiment.
    5. Print the results to show which gene, and therefore which location, has the
       highest frequency.
    """

    # Step 1: Define gene map based on time of entry (in arbitrary minutes)
    # The gene order is thr-azi-gal, so thr enters first.
    gene_entry_times = {
        'thr': 10,
        'azi': 25,
        'gal': 40
    }
    
    # The answer choices map to these locations/markers:
    # A. Immediately after thr+ -> thr locus
    # C. Between azy and gal -> could refer to azi or some marker in between
    # E. Adjacent to gal -> gal locus
    location_map = {
        'A. Immediately after thr+': 'thr',
        'C. Between azy and gal': 'azi', # Representing the next major marker in this region
        'E. Adjacent to gal': 'gal'
    }


    # Step 2: Set up simulation parameters
    num_simulations = 100000
    # Mating duration is random, up to a maximum time.
    max_mating_duration = 50 

    # Step 3: Initialize counters and run the simulation
    transfer_counts = {
        'thr': 0,
        'azi': 0,
        'gal': 0
    }

    for _ in range(num_simulations):
        # Generate a random duration for this mating event
        mating_duration = random.uniform(0, max_mating_duration)

        # Check which genes were transferred based on the duration
        if mating_duration >= gene_entry_times['thr']:
            transfer_counts['thr'] += 1
        if mating_duration >= gene_entry_times['azi']:
            transfer_counts['azi'] += 1
        if mating_duration >= gene_entry_times['gal']:
            transfer_counts['gal'] += 1

    # Step 4 & 5: Calculate and print the frequencies
    print("--- Simulating Frequencies of Recombinants in E. coli ---")
    print(f"Based on gene order thr-azi-gal and {num_simulations} simulated matings.\n")

    for choice, gene in location_map.items():
        count = transfer_counts[gene]
        frequency = count / num_simulations
        
        print(f"Location: \"{choice}\" (Marker: {gene})")
        # Final equation output format requested by the user
        print(f"   - Recombinant Frequency = {count} (transfers) / {num_simulations} (matings) = {frequency:.4f}")

    print("\n--- Conclusion ---")
    print("The highest frequency of recombinants occurs for the gene transferred earliest ('thr').")
    print("This corresponds to the location 'Immediately after thr+'.")


solve_recombination_frequency()
<<<A>>>
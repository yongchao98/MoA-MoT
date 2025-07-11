import random

def solve_interrupted_mating():
    """
    Simulates an Hfr interrupted mating experiment to find the location
    of the highest recombinant frequency.
    """
    # 1. Define the gene order and their time of entry in minutes.
    # The order is thr -> azi -> gal, so thr has the earliest entry time.
    gene_entry_times = {'thr': 2.0, 'azy': 8.0, 'gal': 15.0}
    gene_order = ['thr', 'azy', 'gal']

    # 2. Simulation parameters
    num_simulations = 100000
    # Average time until a random interruption occurs (e.g., spontaneous separation)
    # We use an exponential distribution to model random events over time.
    # The rate (lambda) is 1 / mean_time.
    mean_mating_time = 10.0  # minutes

    # 3. Run the simulation
    # Dictionary to count how many times each gene is successfully transferred
    recombinant_counts = {'thr': 0, 'azy': 0, 'gal': 0}

    for _ in range(num_simulations):
        # Determine the duration of this specific mating by simulating a random interruption
        mating_duration = random.expovariate(1.0 / mean_mating_time)

        # Check which genes were transferred within the simulated time
        if mating_duration >= gene_entry_times['thr']:
            recombinant_counts['thr'] += 1
        if mating_duration >= gene_entry_times['azy']:
            recombinant_counts['azy'] += 1
        if mating_duration >= gene_entry_times['gal']:
            recombinant_counts['gal'] += 1

    # 4. Analyze and print the results
    print("--- Interrupted Mating Simulation ---")
    print(f"Gene transfer order: {' -> '.join(gene_order)}")
    print(f"Entry times (minutes): {gene_entry_times}\n")
    print(f"Simulating {num_simulations} matings with random interruptions...")
    print("--------------------------------------")
    print("Resulting recombinant frequencies:")

    # This loop demonstrates the "final equation" by calculating and showing each frequency.
    # Final Equation: Frequency = (Number of recombinants for gene) / (Total matings)
    for gene in gene_order:
        frequency = recombinant_counts[gene] / num_simulations
        print(f"  - Frequency of {gene}+ = {recombinant_counts[gene]} / {num_simulations} = {frequency:.4f}")

    # 5. Conclusion based on results
    # Find the gene with the highest frequency
    highest_freq_gene = max(recombinant_counts, key=recombinant_counts.get)
    
    print("\n--- Conclusion ---")
    print(f"The simulation shows that '{highest_freq_gene}' has the highest frequency of recombination.")
    print("This is because it is the first gene to be transferred and has the highest probability of entering the recipient cell before mating is interrupted.")
    print("The genetic location of this highest frequency is therefore at the start of the transferred map segment, which is best described as 'Immediately after thr+'.")

solve_interrupted_mating()
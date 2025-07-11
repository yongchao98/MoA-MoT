import math

def simulate_interrupted_mating():
    """
    Simulates an interrupted mating experiment to find the location of the highest
    recombinant frequency.
    """
    # --- Setup ---
    # Gene order on the chromosome
    gene_order = ['thr', 'azi', 'gal']

    # Time of entry (in minutes) for each genetic marker.
    # The first marker, thr, enters early. The subsequent markers enter later.
    entry_times = {
        'thr': 2,
        'azi': 10,
        'gal': 18
    }

    # --- Simulation Parameters ---
    # Let's assume the probability that the mating bridge remains intact for one
    # minute is 0.9 (i.e., a 10% chance of breaking per minute).
    prob_persist_per_minute = 0.9

    # Initial population of mating pairs
    initial_population = 1000

    print("Simulating Hfr Interrupted Mating Experiment")
    print("-" * 50)
    print(f"Gene order: {' -> '.join(gene_order)}")
    print(f"Assumed probability of mating bridge persisting per minute: {prob_persist_per_minute}\n")

    # --- Calculation ---
    print("Calculating expected recombinant frequencies...")
    print("-" * 50)

    recombinant_counts = {}
    for gene in gene_order:
        time = entry_times[gene]
        # The probability of the bridge persisting long enough to transfer a gene is
        # (prob_persist_per_minute) ^ time.
        prob_transfer = prob_persist_per_minute ** time
        
        # The number of recombinants is the initial population * probability of transfer.
        # We assume every successful transfer leads to a recombinant for simplicity.
        num_recombinants = initial_population * prob_transfer
        recombinant_counts[gene] = round(num_recombinants)

    # --- Results and Conclusion ---
    for gene in gene_order:
        count = recombinant_counts[gene]
        time = entry_times[gene]
        frequency = count / initial_population * 100
        print(f"Gene: {gene:<5} | Time of Entry: {time:>2} min | Expected Recombinants: {count:>4} | Frequency: {frequency:.1f}%")

    print("-" * 50)
    print("\nConclusion:")
    print("The simulation shows a gradient of transfer. The earliest gene, 'thr',")
    print("is transferred to the most recipient cells, resulting in the highest")
    print("frequency of recombinants.")
    print("\nTherefore, the genetic location with the highest frequency of recombinants is")
    print("at the 'thr' locus. Among the given choices, 'Immediately after thr+' is the")
    print("best description of this location.")


if __name__ == "__main__":
    simulate_interrupted_mating()

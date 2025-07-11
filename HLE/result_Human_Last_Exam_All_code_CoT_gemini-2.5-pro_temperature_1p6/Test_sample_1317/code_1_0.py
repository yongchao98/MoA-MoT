import random

def simulate_interrupted_mating():
    """
    Simulates Hfr interrupted mating to find the frequency of recombinants.
    """
    # Gene order and their time of entry in minutes (arbitrary but sequential)
    genes = {
        'thr': 8,
        'azi': 17,
        'gal': 25
    }

    # Number of mating pairs to simulate
    num_simulations = 10000

    # Counters for each recombinant type
    recombinant_counts = {gene: 0 for gene in genes}

    # Average time until the mating is naturally interrupted (in minutes)
    # We use this as the rate for an exponential distribution of interruption times.
    # A lambda of 1/20 means an average interruption time of 20 minutes.
    interruption_rate = 1 / 20.0

    print("Starting simulation of Hfr interrupted mating...")
    print(f"Gene order and transfer times: {genes}")
    print("-" * 30)

    for _ in range(num_simulations):
        # Generate a random time when the mating is interrupted.
        # random.expovariate is a good model for random breakage events.
        interruption_time = random.expovariate(interruption_rate)

        # Check which genes were successfully transferred before the interruption.
        # A gene is transferred only if its transfer time is less than the interruption time.
        if interruption_time > genes['thr']:
            recombinant_counts['thr'] += 1
        if interruption_time > genes['azi']:
            recombinant_counts['azi'] += 1
        if interruption_time > genes['gal']:
            recombinant_counts['gal'] += 1
            
    print("Simulation Results:")
    print(f"Total simulated matings: {num_simulations}")
    print("\nFrequency of recombinants for each gene:")
    for gene, count in recombinant_counts.items():
        frequency = (count / num_simulations) * 100
        print(f"- {gene}:\t{count} successful transfers ({frequency:.2f}%)")

    print("\nConclusion:")
    print("The simulation shows that the frequency of recombinants is highest for the gene that is transferred first ('thr') and decreases for subsequent genes.")
    print("Therefore, the highest frequency of recombinants should be expected at the location of the earliest transferred marker.")
    print("Based on the provided choices, the location is immediately after thr+ is transferred.")

if __name__ == '__main__':
    simulate_interrupted_mating()
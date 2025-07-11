import random

def simulate_interrupted_mating():
    """
    Simulates an Hfr interrupted mating experiment to determine which gene
    has the highest frequency of recombination.
    """
    # Define genes and their time of entry in minutes. The order is what matters.
    # thr enters first, so it has the earliest time.
    genes = {
        'thr': {'entry_time': 8, 'recombinants': 0},
        'azi': {'entry_time': 17, 'recombinants': 0},
        'gal': {'entry_time': 25, 'recombinants': 0}
    }

    # Number of simulated mating experiments to run
    num_simulations = 100000

    # In a real experiment, mating is interrupted at various defined times.
    # Here, we simulate this by choosing a random interruption time for each event.
    # We set a max time long enough for all genes to potentially be transferred.
    max_mating_time = 40

    # --- Simulation Start ---
    for _ in range(num_simulations):
        # For each simulation, interrupt mating at a random time.
        interruption_time = random.uniform(0, max_mating_time)

        # Check which genes have been transferred. A gene is considered transferred
        # if the mating duration is longer than or equal to its entry time.
        # The transferred gene can then form a recombinant.
        if interruption_time >= genes['thr']['entry_time']:
            genes['thr']['recombinants'] += 1
        
        if interruption_time >= genes['azi']['entry_time']:
            genes['azi']['recombinants'] += 1
            
        if interruption_time >= genes['gal']['entry_time']:
            genes['gal']['recombinants'] += 1
            
    # --- Simulation End ---

    # --- Analysis and Output ---
    print("--- Interrupted Mating Simulation Results ---")
    print(f"Gene transfer order is known: thr -> azi -> gal")
    print(f"Based on {num_simulations} simulated mating events:\n")
    
    # Calculate and print the frequency for each gene
    for gene, data in genes.items():
        frequency = (data['recombinants'] / num_simulations) * 100
        print(f"Frequency of '{gene}+' recombinants: {data['recombinants']}/{num_simulations} = {frequency:.2f}%")

    print("\n--- Conclusion ---")
    print("The simulation shows that the frequency of recombinants is highest for the 'thr' gene.")
    print("This is because 'thr' is the first gene to be transferred from the donor to the recipient cell.")
    print("Therefore, the genetic location with the highest frequency of recombinants is the locus for the earliest marker, 'thr'.")

# Run the simulation
simulate_interrupted_mating()
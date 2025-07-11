import math

def simulate_recombination_frequency():
    """
    Simulates recombination frequency in an interrupted mating experiment.
    
    The frequency is modeled as a decreasing function of the time of entry,
    reflecting that later genes are less likely to be transferred.
    """
    
    # Gene order and relative time of entry (in minutes)
    genes = {
        'thr': 5,
        'azi': 10,
        'gal': 18
    }
    
    print("Simulating recombination frequency based on time of entry.")
    print("Gene transfer order: thr -> azi -> gal\n")
    print(f"{'Gene':<10} | {'Time of Entry (min)':<20} | {'Relative Frequency (%)':<20}")
    print("-" * 60)
    
    # A simple model: Frequency decreases exponentially with time
    # Frequency = Initial_Frequency * decay_rate ^ time
    # Let's use 100 as initial frequency and 0.9 as decay rate per minute.
    initial_frequency = 100.0
    decay_rate = 0.9
    
    max_freq = 0
    max_freq_gene = ""
    
    for gene, time in genes.items():
        # Calculate frequency
        frequency = initial_frequency * (decay_rate ** time)
        
        # Track the gene with the highest frequency
        if frequency > max_freq:
            max_freq = frequency
            max_freq_gene = gene
            
        print(f"{gene:<10} | {time:<20} | {frequency:<20.2f}")

    print("-" * 60)
    print(f"\nConclusion:")
    print(f"The highest frequency of recombinants is observed for the gene '{max_freq_gene}',")
    print("which is the first gene to be transferred from the donor to the recipient.")
    print("This corresponds to the genetic location immediately after the origin of transfer, starting with thr+.")


simulate_recombination_frequency()
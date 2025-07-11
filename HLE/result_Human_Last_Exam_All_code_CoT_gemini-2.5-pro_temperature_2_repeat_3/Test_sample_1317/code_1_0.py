import sys

def model_recombinant_frequency():
    """
    Models and explains the frequency of recombinants in an Hfr cross.
    """

    # Gene order as determined by the time of entry.
    gene_order = ['thr', 'azi', 'gal']

    # We assume the first gene to be transferred has the highest probability of forming a recombinant.
    # Let's set an initial high frequency.
    # The actual value is arbitrary; its relative size is what matters.
    initial_frequency = 95.0 # Representing 95%

    # As we move along the chromosome, the probability of spontaneous interruption of mating increases.
    # This leads to a lower transfer frequency for genes further from the origin.
    # We'll use an attenuation factor to model this drop.
    attenuation_factor = 0.4 # Each subsequent gene has 40% the frequency of the previous one.

    print("--- Interrupted Mating Recombinant Frequency Model ---")
    print(f"Gene transfer order: {' -> '.join(gene_order)}")
    print("The frequency of recombinants is highest for the gene transferred first")
    print("and decreases for genes transferred later.\n")

    current_frequency = initial_frequency

    # Calculate and print the expected frequency for each gene
    for i, gene in enumerate(gene_order):
        # We need at least one crossover before the gene and one after for stable recombination.
        # The probability is dominated by whether the gene is transferred at all.
        print(f"Expected recombinant frequency at '{gene}' locus: {current_frequency:.1f}%")

        # Apply the attenuation factor for the next gene
        current_frequency *= attenuation_factor
        
    print("\n--- Conclusion ---")
    print("The model demonstrates that the 'thr' gene, being the first to transfer,")
    print("is expected to have the highest frequency of recombinants.")
    sys.stdout.flush()

if __name__ == "__main__":
    model_recombinant_frequency()
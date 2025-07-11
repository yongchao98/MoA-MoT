import sys

def solve_recombination_frequency():
    """
    Models and explains recombination frequency in an interrupted mating experiment.
    """
    # The gene order of transfer is given as thr-azi-gal.
    genes = ['thr', 'azi', 'gal']

    # In conjugation, the frequency of recombination is highest for the earliest
    # transferred gene and decreases for later genes. We can model this with
    # some hypothetical numbers.
    # Let's assume the base frequency for the first gene is 95%.
    base_frequency = 95
    # Let's assume the frequency drops by 30% for each subsequent gene.
    frequency_decay = 30

    print("--- Analysis of Recombination Frequency ---")
    print(f"Gene transfer order: {genes[0]} -> {genes[1]} -> {genes[2]}")
    print("\nThe principle is that recombination frequency is highest for the first gene transferred and decreases for subsequent genes.")
    print("We can model this with a simple calculation:")
    print(f"Frequency(gene) = Base Frequency - (Position * Decay)")
    print(f"Where Base Frequency = {base_frequency}% and Decay = {frequency_decay}% per gene.")
    print("\nCalculated hypothetical frequencies:")

    recombination_frequencies = {}
    highest_freq = -1
    highest_freq_gene = None

    for i, gene in enumerate(genes):
        # Calculate the frequency for the current gene
        # The equation is: frequency = base_frequency - (i * frequency_decay)
        # We'll show each number in the calculation.
        current_frequency = base_frequency - (i * frequency_decay)
        recombination_frequencies[gene] = current_frequency
        
        print(f"For '{gene}' (position {i}):")
        print(f"  Frequency = {base_frequency} - ({i} * {frequency_decay}) = {current_frequency}%")
        
        if current_frequency > highest_freq:
            highest_freq = current_frequency
            highest_freq_gene = gene
    
    print("\n--- Conclusion ---")
    print(f"The model shows the highest frequency of recombinants ({highest_freq}%) is for the '{highest_freq_gene}' gene.")
    print("This is because it is the first marker to be transferred from the donor to the recipient cell.")
    print("Therefore, the genetic location with the highest frequency of recombinants is at the 'thr' locus.")
    print("This corresponds to the answer choice 'Immediately after thr+'.")


# Execute the function
solve_recombination_frequency()

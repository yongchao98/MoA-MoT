import sys

def solve_recombination_frequency():
    """
    Models and explains recombination frequency in an E. coli interrupted mating experiment.
    """
    # The gene order on the chromosome is given as thr-azi-gal.
    gene_order = ['thr', 'azi', 'gal']

    # The time of entry determines recombination frequency. The earlier a gene is
    # transferred, the higher its frequency of recombination. Let's assign
    # illustrative time-of-entry values (in minutes).
    time_of_entry = {
        'thr': 8,   # Transferred first
        'azi': 15,  # Transferred second
        'gal': 25   # Transferred third
    }

    print("--- E. coli Interrupted Mating Model ---")
    print(f"Gene order on chromosome: {gene_order[0]}-{gene_order[1]}-{gene_order[2]}")
    print("Transfer starts at 'thr' and proceeds towards 'gal'.\n")

    print("Principle: Recombination frequency is highest for genes transferred earliest.")
    print("We can model this with the equation: Relative Frequency = Constant / Time of Entry")
    print("Let's use a constant of 1000 for this illustration.\n")

    print("Calculating relative recombination frequencies:")
    
    highest_freq = -1
    highest_freq_gene = ""

    # Calculate and display the frequency for each gene
    for gene in gene_order:
        time = time_of_entry[gene]
        # The equation for this model:
        relative_frequency = 1000 / time
        
        # We must output each number in the final equation as requested.
        # Although there is no single final equation, we show the calculation for each gene.
        print(f"For gene '{gene}': Relative Frequency = 1000 / {time} = {relative_frequency:.2f}")

        if relative_frequency > highest_freq:
            highest_freq = relative_frequency
            highest_freq_gene = gene

    print("\n--- Conclusion ---")
    print(f"The gene '{highest_freq_gene}' is transferred earliest and therefore has the highest recombination frequency.")
    print("This means the genetic location with the highest frequency of recombinants is the one closest to the origin of transfer.")
    print("Based on the data, this location is 'thr' and the region immediately following it.")

# Execute the function to print the explanation.
solve_recombination_frequency()
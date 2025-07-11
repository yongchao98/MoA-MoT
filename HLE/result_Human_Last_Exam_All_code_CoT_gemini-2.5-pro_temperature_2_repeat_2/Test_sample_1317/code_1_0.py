import sys

def solve_recombination_frequency():
    """
    This script explains and demonstrates why the highest frequency of recombinants
    in an interrupted mating experiment occurs for the gene transferred first.
    """

    # 1. Define the genes and their established order of entry.
    # We assign arbitrary time units to represent the time of entry. A lower number means earlier entry.
    gene_entry_times = {
        'thr': 10,
        'azy': 25,
        'gal': 40
    }
    gene_order = ['thr', 'azy', 'gal']

    # 2. Explain the core principle.
    print("--- Principle ---")
    print("In interrupted mating, genes are transferred sequentially from an Hfr cell to an F- cell.")
    print("The frequency of recombination for a gene is highest for the genes that are transferred earliest.")
    print("This is because they have the most time to enter the recipient cell before mating is interrupted.\n")
    print(f"Given gene order: {' -> '.join(gene_order)}\n")


    # 3. Calculate and display a hypothetical frequency for each gene.
    # We model frequency as being inversely proportional to the entry time.
    # Equation: Relative Frequency = Constant / Entry Time
    # We will use a constant of 1000 for this demonstration.
    print("--- Calculation ---")
    print("We can model this with the equation: Relative Frequency = 1000 / Entry Time")

    highest_frequency = -1.0
    location_of_highest_frequency = None

    for gene in gene_order:
        time = gene_entry_times[gene]
        # Calculate the frequency using our model.
        frequency = 1000 / time

        # Keep track of the highest frequency found so far.
        if frequency > highest_frequency:
            highest_frequency = frequency
            location_of_highest_frequency = gene

        # Output the equation with the numbers for each gene, as requested.
        print(f"For gene '{gene}': 1000 / {time} = {frequency:.2f}")

    # 4. State the conclusion based on the results.
    print("\n--- Conclusion ---")
    print(f"The first gene transferred, '{location_of_highest_frequency}', shows the highest relative recombination frequency ({highest_frequency:.2f}).")
    print("Therefore, in the given experiment, the highest frequency of recombinants should be observed at the location of the 'thr+' marker, which is the first one to enter the recipient cell.")

# Execute the function to solve the problem.
solve_recombination_frequency()
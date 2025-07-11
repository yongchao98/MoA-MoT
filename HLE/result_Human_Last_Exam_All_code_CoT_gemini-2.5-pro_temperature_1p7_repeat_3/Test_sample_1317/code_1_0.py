import math

def calculate_recombination_frequency():
    """
    This function models the recombination frequency in an interrupted mating experiment.
    The principle is that the frequency of recombination is highest for genes that are
    transferred earliest into the recipient cell.
    """

    # Step 1: Define the gene order and their hypothetical entry times in minutes.
    # 'thr' is the first gene transferred, as stated in the problem.
    genes_and_times = {
        'thr': 9,
        'azi': 18,
        'gal': 26
    }
    
    # Step 2: Define a decay constant 'k' for the model.
    # The probability of the mating bridge remaining intact until time 't' can be
    # modeled as P(t) = exp(-k * t). This probability is proportional to the
    # recombination frequency. A larger time 't' results in a lower probability.
    k = 0.05

    print("Modeling Recombination Frequency based on Time of Gene Entry\n")
    print(f"Gene order: thr -> azi -> gal")
    print(f"Model for transfer success probability (P) at time (t): P = exp(-k * t), where k = {k}\n")
    
    recombination_frequencies = {}

    # Step 3: Calculate the relative recombination frequency for each gene.
    # We will print out the calculation for each gene.
    for gene, time in genes_and_times.items():
        # The equation for the relative frequency
        frequency = math.exp(-k * time)
        recombination_frequencies[gene] = frequency
        
        # Outputting each number in the final equation as requested
        print(f"For gene '{gene}':")
        print(f"  Entry Time (t) = {time} minutes")
        print(f"  Relative Frequency = exp(-{k} * {time}) = {frequency:.4f}")

    # Step 4: Determine and state the conclusion.
    max_freq_gene = max(recombination_frequencies, key=recombination_frequencies.get)

    print("\n---")
    print(f"Conclusion: The model shows that the gene '{max_freq_gene}' has the highest relative recombination frequency.")
    print("This is because it is the first gene to enter the recipient cell, giving it the most time and highest probability to recombine before the mating is interrupted.")
    print("Therefore, the highest frequency of recombinants should be expected at the location of the thr+ gene.")

# Run the simulation
calculate_recombination_frequency()
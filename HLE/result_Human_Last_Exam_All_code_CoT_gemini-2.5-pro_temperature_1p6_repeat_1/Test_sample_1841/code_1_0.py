def calculate_simpsons_diversity():
    """
    Calculates and explains Simpson's Diversity Index for a sample with one species.
    """
    # Sample data: A single species of bat is found, with 75 individuals.
    species_counts = [75]
    
    # N is the total number of organisms of all species.
    N = sum(species_counts)
    
    # sum_n_minus_1 is the sum of n(n-1) for each species.
    sum_n_minus_1 = 0
    for n in species_counts:
        sum_n_minus_1 += n * (n - 1)
        
    # N_N_minus_1 is N(N-1)
    # Check to prevent division by zero if N is 0 or 1.
    if N < 2:
        print("Simpson's Index is not defined for samples with fewer than 2 individuals.")
        return
        
    N_N_minus_1 = N * (N - 1)
    
    # Simpson's Diversity Index D
    D = 1 - (sum_n_minus_1 / N_N_minus_1)
    
    print("This script demonstrates that a Simpson's Diversity Index of 0 is mathematically possible.")
    print("It occurs when only one species is present in the sample.")
    print("\n--- Calculation ---")
    print(f"Sampled Species Counts (n): {species_counts}")
    print(f"Total Individuals (N): {N}")
    print(f"Sum of n(n-1): {sum_n_minus_1}")
    print(f"N(N-1): {N_N_minus_1}")
    print("\n--- Final Equation ---")
    print(f"D = 1 - ( {sum_n_minus_1} / {N_N_minus_1} )")
    print(f"D = 1 - {sum_n_minus_1 / N_N_minus_1}")
    print(f"D = {D}")
    print("\nThis result is mathematically valid. However, it is not ecologically valid because the problem states several bat species are known to live on the island, meaning the sample was not representative of the true diversity.")

calculate_simpsons_diversity()
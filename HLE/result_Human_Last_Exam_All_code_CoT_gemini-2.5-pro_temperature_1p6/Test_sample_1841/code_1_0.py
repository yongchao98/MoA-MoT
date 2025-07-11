import math

def calculate_simpsons_diversity():
    """
    Demonstrates the calculation of Simpson's Diversity Index for a sample with only one species.
    D = 1 - (sum(n*(n-1))) / (N*(N-1))
    """

    print("Analyzing the Simpson's Diversity Index result of 0.")
    print("----------------------------------------------------\n")
    print("Let's assume the student's sample consists of 25 individuals from a single bat species.")

    # n_i = number of individuals of species i
    # In this case, we only have one species.
    species_counts = [25]
    species_names = ["Bat Species A"]

    # N = total number of individuals
    N = sum(species_counts)

    # sum_n_minus_1 = sum of n*(n-1) for each species
    sum_n_minus_1 = 0
    for n in species_counts:
        sum_n_minus_1 += n * (n - 1)

    # N_minus_1 = N * (N - 1)
    # Check for N > 1 to avoid division by zero, although D is undefined for N<=1
    if N <= 1:
        print("Simpson's Index is not defined for a total population (N) of 1 or less.")
        return

    N_minus_1 = N * (N - 1)
    
    # Calculate the Diversity Index (D)
    # Handle the case where the denominator could be zero
    if N_minus_1 == 0:
        # This occurs if N=0 or N=1, meaning no diversity by definition.
        simpsons_index_D = 0
    else:
        dominance_index_lambda = sum_n_minus_1 / N_minus_1
        simpsons_index_D = 1 - dominance_index_lambda
    
    print("\n--- Calculation Steps ---")
    print(f"Species found: {species_names[0]}")
    print(f"Number of individuals of this species (n): {species_counts[0]}")
    print(f"Total number of individuals in the sample (N): {N}")
    print("\nFormula: D = 1 - (Σn(n-1)) / (N(N-1))\n")
    
    print(f"Calculating the numerator part, Σn(n-1):")
    # This loop shows the calculation for each term, even though there's only one.
    sum_n_parts_str = []
    for i, n in enumerate(species_counts):
        calc = n * (n - 1)
        print(f"  For {species_names[i]}: n(n-1) = {n}({n}-1) = {calc}")
        sum_n_parts_str.append(str(calc))
    print(f"Σn(n-1) = {' + '.join(sum_n_parts_str)} = {sum_n_minus_1}\n")

    print(f"Calculating the denominator part, N(N-1):")
    print(f"N(N-1) = {N}({N}-1) = {N_minus_1}\n")

    print(f"Putting it all together:")
    print(f"D = 1 - {sum_n_minus_1} / {N_minus_1}")
    print(f"D = 1 - {sum_n_minus_1 / N_minus_1}")
    print(f"D = {simpsons_index_D}\n")
    
    print("--- Conclusion ---")
    print("The calculation shows that a Simpson's Diversity Index of 0 is mathematically valid and occurs when only one species is found.")
    print("However, finding only one species on an island known to have several makes the result ecologically questionable or unrepresentative of the true ecosystem.")
    print("Therefore, the index value is 'Mathematically valid, but not ecologically valid.'")

calculate_simpsons_diversity()
def calculate_simpsons_index(species_counts):
    """
    Calculates Simpson's Diversity Index (D) from a dictionary of species counts.
    The formula is D = 1 - (Σn(n-1)) / (N(N-1))
    
    Args:
        species_counts (dict): A dictionary where keys are species names
                               and values are the number of individuals.
    """
    
    # Filter out species with zero counts to represent what was actually found
    found_species = {species: n for species, n in species_counts.items() if n > 0}
    
    # If no species or only one species is found, the logic proceeds.
    # N is the total number of organisms of all species found.
    N = sum(found_species.values())
    
    # The sum of n(n-1) for each species.
    # 'n' is the number of organisms of a particular species.
    sum_n_n_minus_1 = sum(n * (n - 1) for n in found_species.values())
    
    # The total number of organisms multiplied by (total - 1).
    N_N_minus_1 = N * (N - 1)
    
    # If N(N-1) is zero (i.e., total individuals N is 0 or 1), 
    # diversity is conventionally 0.
    if N_N_minus_1 == 0:
        D = 0.0
        print("The sample contains 1 or 0 individuals, so diversity is 0 by definition.")
        return

    # Calculate Simpson's Index D
    D = 1 - (sum_n_n_minus_1 / N_N_minus_1)
    
    print("Demonstrating the calculation for a sample with zero diversity:")
    print(f"Sampled Counts: {species_counts}")
    print("-" * 30)
    print(f"Total individuals found (N): {N}")
    print(f"Sum of n(n-1) for each species [Σn(n-1)]: {sum_n_n_minus_1}")
    print(f"Total individuals calculation [N(N-1)]: {N_N_minus_1}")
    print("\nFinal Equation:")
    print(f"D = 1 - ( {sum_n_n_minus_1} / {N_N_minus_1} )")
    print(f"D = 1 - {sum_n_n_minus_1 / N_N_minus_1}")
    print(f"D = {D}")
    print("\nThis shows that D=0 is a mathematically valid result for a sample containing only one species.")


# --- Main execution ---
# Let's simulate the student's sample. Assume they went out and only found
# individuals of a single bat species, for example, 30 individuals of 'Bat Species A'.
student_sample = {
    'Bat Species A': 30,
    'Bat Species B': 0,
    'Bat Species C': 0
}

calculate_simpsons_index(student_sample)
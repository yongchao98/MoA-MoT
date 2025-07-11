import math

def calculate_and_print_shannon():
    """
    This script calculates the Shannon's diversity for three insect survey sites
    after accounting for the true species identity based on visual analysis.
    """
    
    # Step 1: Define the relationship between morphospecies and true species.
    # Based on visual analysis, A, B, C, D, and E are all color morphs of the same species.
    # F is a different species.
    morpho_to_true = {
        'A': 'Species 1 (Harmonia axyridis)',
        'B': 'Species 1 (Harmonia axyridis)',
        'C': 'Species 1 (Harmonia axyridis)',
        'D': 'Species 1 (Harmonia axyridis)',
        'E': 'Species 1 (Harmonia axyridis)',
        'F': 'Species 2'
    }

    # Step 2: Define the raw survey data for each site.
    site_data = {
        "Site 1": {'A': 1, 'B': 1, 'C': 2},
        "Site 2": {'A': 10, 'D': 12, 'E': 1},
        "Site 3": {'B': 2, 'E': 2, 'F': 4}
    }

    print("Based on the images, morphospecies A, B, C, D, and E are the same species.")
    print("Recalculating species counts and computing Shannon's Diversity (H') for each site.")
    print("Formula: H' = -Σ(pi * ln(pi)), where pi = (count of a species / total individuals)\n")

    shannon_results = []

    # Step 3 & 4: Process each site
    for site_name, morpho_counts in site_data.items():
        # Consolidate counts based on true species identity
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = morpho_to_true[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        # Calculate total individuals (N)
        N = sum(true_species_counts.values())
        
        shannon_index = 0
        equation_parts = []
        
        if N == 0 or len(true_species_counts) <= 1:
            # If only one species, diversity is 0 (ln(1) = 0)
            shannon_index = 0
            if N > 0:
                p = N / N
                equation_parts.append(f"{p:.1f} * ln({p:.1f})")
        else:
            # Calculate Shannon index for multiple species
            for count in true_species_counts.values():
                p_i = count / N
                shannon_index += p_i * math.log(p_i)
                equation_parts.append(f"({count}/{N} * ln({count}/{N}))")

        # Final diversity value is the negative sum
        final_shannon = -shannon_index
        shannon_results.append(final_shannon)
        
        # Print the calculation steps for the current site
        print(f"{site_name}:")
        print(f"  True species counts: {true_species_counts}")
        equation_str = " + ".join(equation_parts)
        print(f"  H' = - [ {equation_str} ] = {final_shannon:.2f}\n")

    # Step 5: Format and print the final consolidated answer
    final_answer_string = ", ".join([f"{val:.2f}" for val in shannon_results])
    print("Final Answer:")
    print(final_answer_string)

calculate_and_print_shannon()
<<<0.00, 0.00, 0.69>>>
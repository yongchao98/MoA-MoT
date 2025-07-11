import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites after
    correcting for species identity based on polymorphism.
    """
    # Step 1: Define the initial abundance data for each morphospecies at each site.
    sites_morpho_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    # Step 2: Define the mapping from morphospecies to true species.
    # A, B, C, D, E are all forms of Harmonia axyridis (Species 1).
    # F is a different species (Species 2).
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis', 'F': 'Species 2'
    }

    results = []
    print("Calculating insect diversity for each site...\n")

    # Step 3 & 4: Process each site to calculate diversity.
    for site_num, morpho_counts in sites_morpho_data.items():
        print(f"--- Site {site_num} ---")
        
        # Consolidate counts based on true species identity
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        print(f"Corrected Abundances: {true_species_counts}")

        total_individuals = sum(true_species_counts.values())
        
        shannon_index = 0.0
        equation_parts = []

        # Handle cases with 0 or 1 species, where diversity is 0.
        if len(true_species_counts) <= 1:
            shannon_index = 0.0
            if total_individuals > 0:
                p = 1.0
                equation_parts.append(f"({p:.2f} * ln({p:.2f}))")
        else:
            for species, count in true_species_counts.items():
                if count > 0:
                    proportion = count / total_individuals
                    shannon_index -= proportion * math.log(proportion)
                    equation_parts.append(f"({proportion:.2f} * ln({proportion:.2f}))")

        # Print the equation used for the calculation
        equation_str = " + ".join(equation_parts)
        print(f"Equation: H = -[ {equation_str} ]")
        print(f"Shannon's Diversity (H): {shannon_index:.2f}\n")
        results.append(shannon_index)

    # Final formatted output
    final_answer_str = ", ".join([f"{res:.2f}" for res in results])
    print("The insect diversity of each site is:")
    print(final_answer_str)
    print(f"<<<{final_answer_str}>>>")

# Run the calculation
calculate_shannon_diversity()
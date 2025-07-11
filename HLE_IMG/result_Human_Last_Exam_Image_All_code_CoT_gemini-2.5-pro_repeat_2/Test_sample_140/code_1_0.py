import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity index for three insect survey sites
    after correcting for true species identities.
    """
    
    # Abundance of each morphospecies at each site
    site_morpho_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    # Mapping from morphospecies to true species
    # Morphospecies A, B, C, and E are all Harmonia axyridis (True Species 1)
    # Morphospecies D is True Species 2
    # Morphospecies F is True Species 3
    species_map = {
        'A': 'True Species 1 (H. axyridis)',
        'B': 'True Species 1 (H. axyridis)',
        'C': 'True Species 1 (H. axyridis)',
        'D': 'True Species 2',
        'E': 'True Species 1 (H. axyridis)',
        'F': 'True Species 3'
    }

    print("Step 1: Re-calculating species abundance based on true species identities.\n")
    print("Morphospecies A, B, C, and E are the same species (H. axyridis).\n")

    results = []
    
    for site_num, morpho_counts in sorted(site_morpho_data.items()):
        print(f"--- Calculating Diversity for Site {site_num} ---")
        
        # Calculate true species counts for the site
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        print(f"True species counts: {true_species_counts}")
        
        total_individuals = sum(true_species_counts.values())
        print(f"Total individuals: {total_individuals}")
        
        shannon_index = 0.0
        
        equation_parts = []
        for species, count in true_species_counts.items():
            if count > 0:
                pi = count / total_individuals
                shannon_index -= pi * math.log(pi)
                equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
        
        # Print the full equation for clarity
        if not equation_parts:
             print("Shannon's Diversity Equation: H = 0 (No individuals found)")
        else:
             print(f"Shannon's Diversity Equation: H = - [ {' + '.join(equation_parts)} ]")

        results.append(shannon_index)
        print(f"Shannon's Diversity (H) for Site {site_num}: {shannon_index:.2f}\n")

    # Final formatted answer
    final_answer = ", ".join([f"{r:.2f}" for r in results])
    print("--- Final Answer ---")
    print("The insect diversity of each site is:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


calculate_shannon_diversity()
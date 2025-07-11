import math

def solve_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites after correcting
    species identification from morphospecies to true species.
    """
    print("Step 1: Define initial data and taxonomic corrections.")
    # Abundance of each morphospecies at each site.
    site_data_morpho = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    # Based on the images, morphospecies A, B, C, D, and E are all identified as the same
    # polymorphic species, Harmonia axyridis. Morphospecies F is a distinct species.
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis', 'F': 'Distinct Species F'
    }
    print("Taxonomic correction: Morphospecies A, B, C, D, E are one species. F is another.")
    print("-" * 40)

    shannon_values = []
    
    print("Step 2: Calculate true species abundance and Shannon's diversity for each site.")
    # Process each site
    for site_num in sorted(site_data_morpho.keys()):
        print(f"\nProcessing Site {site_num}:")
        morpho_counts = site_data_morpho[site_num]
        
        # Consolidate abundances based on true species
        true_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_counts[true_species] = true_counts.get(true_species, 0) + count
            
        print(f"  - True species abundances: {true_counts}")
        
        abundances = list(true_counts.values())
        total_individuals = sum(abundances)
        
        # Calculate Shannon's Diversity H' = -Σ(pi * ln(pi))
        shannon_h = 0.0
        if total_individuals == 0 or len(abundances) <= 1:
            # Diversity is 0 if there is 0 or 1 species
            shannon_h = 0.0
            print(f"  - Only one species type found. Shannon's diversity (H') is 0.00.")
        else:
            print(f"  - Shannon's Diversity H' = -Σ (pi * ln(pi))")
            equation_parts = []
            term_values = []
            for abundance in abundances:
                proportion = abundance / total_individuals
                term = proportion * math.log(proportion)
                term_values.append(term)
                equation_parts.append(f"({abundance}/{total_individuals} * ln({abundance}/{total_individuals}))")
            
            shannon_h = -sum(term_values)
            
            # Show the equation with numbers
            equation_str = " - (" + " + ".join(equation_parts) + ")"
            print(f"  - H' = {equation_str}")
            print(f"  - H' = - ({' + '.join([f'{t:.3f}' for t in term_values])})")
            print(f"  - H' = - ({sum(term_values):.3f})")
            print(f"  - H' = {shannon_h:.2f}")

        shannon_values.append(shannon_h)
    
    print("-" * 40)
    
    # Format and print the final answer string
    final_answer_str = ", ".join([f"{val:.2f}" for val in shannon_values])
    print("\nFinal Answer:")
    print("The insect diversity of each site is:")
    print(final_answer_str)

solve_diversity()
<<<0.00, 0.00, 0.69>>>
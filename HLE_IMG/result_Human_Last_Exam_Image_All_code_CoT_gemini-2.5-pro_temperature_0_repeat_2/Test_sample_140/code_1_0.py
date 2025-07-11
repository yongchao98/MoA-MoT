import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites after
    re-classifying morphospecies into true species.
    """
    # Step 1: Define the initial abundance data for morphospecies at each site.
    site_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    # Step 2: Define the mapping from morphospecies to true species.
    # Based on visual identification:
    # Species 1: Harmonia axyridis (A, B, C, E)
    # Species 2: Coccinella sp. (D)
    # Species 3: Brachiacantha ursina (F)
    species_map = {
        'A': 'Species 1', 'B': 'Species 1', 'C': 'Species 1',
        'D': 'Species 2', 'E': 'Species 1', 'F': 'Species 3'
    }

    results = []
    print("Calculating Shannon's Diversity for each site based on true species identity.\n")

    # Loop through each site to perform calculations
    for site_num, morpho_counts in sorted(site_data.items()):
        
        # Step 3: Re-calculate abundances based on true species.
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        abundances = list(true_species_counts.values())
        total_individuals = sum(abundances)
        
        # Step 4: Calculate Shannon's Diversity Index (H).
        shannon_index = 0
        equation_parts = []
        
        if total_individuals > 0:
            for count in abundances:
                if count > 0:
                    proportion = count / total_individuals
                    shannon_index -= proportion * math.log(proportion)
                    equation_parts.append(f"({count}/{total_individuals} * ln({count}/{total_individuals}))")

        # Print the detailed calculation for the current site
        equation_str = " + ".join(equation_parts)
        print(f"Site {site_num} Abundances: {true_species_counts}")
        print(f"Site {site_num} Shannon's Diversity (H) = -[ {equation_str} ] = {shannon_index:.2f}\n")
        
        results.append(f"{shannon_index:.2f}")

    # Final formatted answer
    final_answer = ", ".join(results)
    print("The insect diversity of each site is:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")

calculate_shannon_diversity()
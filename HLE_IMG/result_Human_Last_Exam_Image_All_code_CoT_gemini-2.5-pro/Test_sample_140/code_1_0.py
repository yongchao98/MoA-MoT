import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites
    after consolidating morphospecies into true species.
    """
    # Step 1: Define initial survey data and species mapping
    site_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    species_map = {
        'A': 'Harmonia axyridis',
        'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis',
        'E': 'Harmonia axyridis',
        'F': 'Tortoise Beetle'
    }

    results = []
    
    print("Calculating Shannon's Diversity for each site:\n")

    # Step 2 & 3: Iterate through sites, recalculate abundances, and calculate diversity
    for site_id in sorted(site_data.keys()):
        morpho_counts = site_data[site_id]
        
        # Recalculate counts based on true species
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count

        # Calculate Shannon's Diversity
        total_individuals = sum(true_species_counts.values())
        shannon_index = 0
        
        equation_parts = []
        if total_individuals > 0:
            for species, count in true_species_counts.items():
                if count > 0:
                    proportion = count / total_individuals
                    shannon_index -= proportion * math.log(proportion)
                    equation_parts.append(f"({count}/{total_individuals} * ln({count}/{total_individuals}))")

        # Print the calculation steps
        equation_str = " + ".join(equation_parts)
        print(f"Site {site_id}:")
        print(f"Species Counts: {true_species_counts}")
        if shannon_index == 0:
             print(f"H' = 0 (Only one species present)")
        else:
             print(f"H' = -[ {equation_str} ]")
        print(f"Result: {shannon_index:.2f}\n")
        
        results.append(f"{shannon_index:.2f}")

    # Final formatted answer
    final_answer = ", ".join(results)
    print("Final Answer:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


calculate_shannon_diversity()
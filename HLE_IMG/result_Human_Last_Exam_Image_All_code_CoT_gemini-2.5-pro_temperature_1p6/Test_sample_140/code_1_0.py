import math

def calculate_shannon_diversity():
    """
    Identifies true species, recalculates abundances, and computes
    Shannon's diversity for three insect survey sites.
    """
    
    # Step 1: Define original abundance data and true species mapping
    site_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    species_map = {
        'A': 'Harmonia axyridis',
        'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis',
        'D': 'Coccinella septempunctata',
        'E': 'Harmonia axyridis',
        'F': 'Not an Insect'
    }

    print("Step 1: Identifying true species and recalculating abundances.\n")
    print("Based on taxonomy:")
    print("- Morphospecies A, B, C, E are identified as one species: Harmonia axyridis.")
    print("- Morphospecies D is identified as a second species: Coccinella septempunctata.")
    print("- Morphospecies F is identified as a spider, which is not an insect, and is excluded.\n")

    diversity_results = []
    
    # Steps 2 & 3: Process each site
    for site_num in sorted(site_data.keys()):
        site_morphospecies = site_data[site_num]
        
        true_species_counts = {}
        # Aggregate counts by true species, ignoring non-insects
        for morpho, count in site_morphospecies.items():
            true_species = species_map[morpho]
            if true_species != 'Not an Insect':
                true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count

        counts_for_calc = list(true_species_counts.values())
        total_individuals = sum(counts_for_calc)

        print(f"--- Calculating Diversity for Site {site_num} ---")
        print(f"Original morphospecies counts: {site_morphospecies}")
        print(f"Adjusted insect species counts: {true_species_counts}")

        shannon_index = 0.0
        # If there are no individuals or only one species, diversity is 0
        if not counts_for_calc or len(counts_for_calc) <= 1:
             print("Since there is only one insect species, the diversity is 0.\n")
        else:
            equation_parts = []
            proportions = []
            for count in counts_for_calc:
                proportion = count / total_individuals
                proportions.append(proportion)
                # Check to avoid math domain error with log(0) if a proportion is 0
                if proportion > 0:
                    shannon_index -= proportion * math.log(proportion)
                    equation_parts.append(f"({count}/{total_individuals} * ln({count}/{total_individuals}))")

            equation = "H = -[" + " + ".join(equation_parts) + "]"
            print(f"Total insect individuals (N) = {total_individuals}")
            print(f"Shannon diversity equation: {equation}")
            print(f"Calculated Shannon's diversity (H) = {shannon_index:.2f}\n")
        
        diversity_results.append(shannon_index)

    # Step 4: Format and print the final answer
    formatted_results = ", ".join([f"{result:.2f}" for result in diversity_results])
    
    print("--------------------------------------------------")
    print(f"The insect diversity of each site is: {formatted_results}")
    return formatted_results

final_answer = calculate_shannon_diversity()
print(f"<<<{final_answer}>>>")

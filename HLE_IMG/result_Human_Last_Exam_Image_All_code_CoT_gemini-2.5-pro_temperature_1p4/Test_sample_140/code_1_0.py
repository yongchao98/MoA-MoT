import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites
    after identifying the true species from morphospecies.
    """
    # Step 1: Define the mapping from morphospecies to true species.
    # Based on the image, A, B, C, D, and E are all morphs of Harmonia axyridis.
    # F is a different species.
    species_map = {
        'A': 'Harmonia axyridis',
        'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis',
        'E': 'Harmonia axyridis',
        'F': 'Species F'
    }

    # Step 2: Define the initial abundance data for each site.
    site_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    results = []

    print("Calculating Shannon's diversity for each site...\n")

    # Step 3: Iterate through each site, aggregate data, and calculate diversity.
    for site_number, abundances in sorted(site_data.items()):
        print(f"--- Site {site_number} ---")
        
        # Aggregate counts based on true species
        true_species_counts = {}
        for morpho, count in abundances.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        print(f"Aggregated counts: {true_species_counts}")
        
        total_individuals = sum(true_species_counts.values())
        print(f"Total individuals (N): {total_individuals}")
        
        shannon_index = 0.0
        
        if total_individuals > 0 and len(true_species_counts) > 1:
            # The calculation is only non-zero if there is more than one species.
            equation_parts = []
            for species, count in true_species_counts.items():
                proportion = count / total_individuals
                shannon_index -= proportion * math.log(proportion)
                equation_parts.append(f"({proportion:.2f} * ln({proportion:.2f}))")
            
            print(f"Calculation: H = -[ {' + '.join(equation_parts)} ]")

        elif len(true_species_counts) == 1:
            # If only one species, diversity is 0.
            p = 1.0
            print(f"Calculation: H = -[ ({p:.2f} * ln({p:.2f})) ] = 0")
        
        else: # No individuals
            print("No individuals found, diversity is 0.")

        print(f"Shannon's Diversity (H) for Site {site_number}: {shannon_index:.2f}\n")
        results.append(f"{shannon_index:.2f}")

    # Step 4: Format and print the final answer string.
    final_answer = ", ".join(results)
    print("--- Final Answer ---")
    print("The insect diversity of each site is:")
    print(final_answer)
    return final_answer

# Execute the function and capture the final answer for the '<<<>>>' format.
final_answer_string = calculate_shannon_diversity()

print(f"\n<<< {final_answer_string} >>>")
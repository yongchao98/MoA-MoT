import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three sites after re-identifying species.
    """
    # Step 1: Define morphospecies to true species mapping
    # Morphospecies A, B, C, D, E are all forms of Harmonia axyridis (True Species 1)
    # Morphospecies F is a different species (True Species 2)
    morpho_to_true_species = {
        'A': 'True Species 1',
        'B': 'True Species 1',
        'C': 'True Species 1',
        'D': 'True Species 1',
        'E': 'True Species 1',
        'F': 'True Species 2'
    }

    # Original abundance data
    site_data = {
        'Site 1': {'A': 1, 'B': 1, 'C': 2},
        'Site 2': {'A': 10, 'D': 12, 'E': 1},
        'Site 3': {'B': 2, 'E': 2, 'F': 4}
    }

    shannon_indices = []
    print("Calculating Shannon's Diversity for each site:\n")

    # Loop through each site to calculate H
    for site_name, morpho_counts in site_data.items():
        # Step 2: Recalculate abundances based on true species
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = morpho_to_true_species[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count

        total_individuals = sum(true_species_counts.values())
        
        print(f"--- {site_name} ---")
        print(f"Abundances of true species: {true_species_counts}")

        # Step 3: Calculate Shannon's diversity index (H)
        shannon_sum = 0
        
        # If there is only one species, diversity is 0
        if len(true_species_counts) < 2:
            shannon_h = 0.0
            print("Calculation: Site has only one species. H = 0\n")
        else:
            equation_parts = []
            for species, count in true_species_counts.items():
                proportion = count / total_individuals
                if proportion > 0:
                    shannon_sum += proportion * math.log(proportion)
                    # Show each number in the equation
                    equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
            
            shannon_h = -shannon_sum
            equation_str = " + ".join(equation_parts)
            print(f"Calculation: H = -[ {equation_str} ]")
            print(f"Result: H = {shannon_h:.2f}\n")
            
        shannon_indices.append(shannon_h)

    # Format final answer
    final_answer = ", ".join([f"{h:.2f}" for h in shannon_indices])
    print("--- Final Answer ---")
    print("The Shannon's diversity of insects at each site is:")
    print(final_answer)

calculate_shannon_diversity()
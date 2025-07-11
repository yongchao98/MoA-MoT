import math

def solve_diversity():
    """
    This function calculates and prints the Shannon's diversity for three insect survey sites
    after accounting for polymorphic species.
    """
    # Step 1: Define initial data from the problem description
    site_data = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }

    # Based on visual analysis, A-E are the same species (Harmonia axyridis), F is different.
    species_mapping = {
        'A': 'Harmonia axyridis',
        'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis',
        'E': 'Harmonia axyridis',
        'F': 'Species F'
    }

    # List to store the final diversity values for each site
    diversities = []

    print("Calculating insect diversity for each site:\n")

    # Step 2 & 3: Process each site
    for site_id in sorted(site_data.keys()):
        morpho_counts = site_data[site_id]
        
        # Aggregate counts by true species
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_mapping[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count

        print(f"--- Site {site_id} ---")
        print(f"True species abundances: {true_species_counts}")

        abundances = list(true_species_counts.values())
        total_individuals = sum(abundances)
        num_species = len(abundances)

        # Calculate Shannon's Diversity (H)
        shannon_index = 0.0
        
        print("Shannon's Diversity Calculation (H = -Σ(pᵢ * ln(pᵢ))):")
        if num_species <= 1:
            shannon_index = 0.0
            print("  Since there is only one species, diversity H = 0.00")
        else:
            equation_parts = []
            for count in abundances:
                proportion = count / total_individuals
                shannon_index -= proportion * math.log(proportion)
                # Format the numbers for the equation string
                p_str = f"{count}/{total_individuals}"
                equation_parts.append(f"({p_str} * ln({p_str}))")
            
            equation = "H = -[" + " + ".join(equation_parts) + "]"
            print(f"  Equation: {equation}")
            print(f"  H = {shannon_index:.2f}")

        diversities.append(f"{shannon_index:.2f}")
        print("-" * 17 + "\n")

    # Final formatted output
    final_answer = ", ".join(diversities)
    print(f"The Shannon's diversity of insects at each site is: {final_answer}")
    return final_answer

# Execute the function to get the answer
final_result = solve_diversity()
# The final answer in the required format
# print(f"<<<{final_result}>>>")
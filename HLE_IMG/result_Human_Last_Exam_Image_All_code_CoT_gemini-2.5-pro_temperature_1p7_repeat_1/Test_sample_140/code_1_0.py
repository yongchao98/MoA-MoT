import math

def calculate_shannon_diversity():
    """
    This function solves the entomology diversity problem by:
    1. Defining the observed morphospecies counts per site.
    2. Consolidating morphospecies into true species based on biology.
    3. Calculating true species abundances for each site.
    4. Calculating and printing the Shannon's diversity index for each site, including the formula.
    5. Printing the final formatted answer.
    """
    # Step 1: Original abundance data of morphospecies at each site
    sites_morpho_data = {
        "Site 1": {'A': 1, 'B': 1, 'C': 2},
        "Site 2": {'A': 10, 'D': 12, 'E': 1},
        "Site 3": {'B': 2, 'E': 2, 'F': 4}
    }

    # Step 2: Define true species based on morphospecies identification
    # Morphospecies A, B, C, D, E are all color morphs of Harmonia axyridis.
    # Morphospecies F is a distinct species (Brachiacantha ursina).
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis',
        'F': 'Brachiacantha ursina'
    }

    # Step 3: Recalculate abundances based on true species for each site
    sites_true_abundances = {}
    for site, morpho_counts in sites_morpho_data.items():
        true_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_counts[true_species] = true_counts.get(true_species, 0) + count
        sites_true_abundances[site] = true_counts

    # Step 4: Calculate and print Shannon's Diversity for each site
    print("Calculating Shannon's Diversity (H) for each site.")
    print("Formula: H = - sum(pi * ln(pi)), where pi = (count of species / total individuals)\n")
    
    diversity_results = []
    # Sort sites for consistent output order
    sorted_sites = sorted(sites_true_abundances.keys())

    for site in sorted_sites:
        abundances = sites_true_abundances[site]
        print(f"--- {site} ---")
        print(f"Species Abundances: {abundances}")
        
        counts = list(abundances.values())
        total_individuals = sum(counts)
        
        if total_individuals == 0 or len(counts) <= 1:
            shannon_index = 0.0
            p1 = 1.0
            n1 = total_individuals
            equation_str = f"H = - [ ({n1}/{total_individuals}) * ln({n1}/{total_individuals}) ]"
            print(equation_str)
            print("Since there is only one species, diversity is 0.")

        else:
            shannon_index = 0
            equation_parts = []
            for count in counts:
                pi = count / total_individuals
                shannon_index -= pi * math.log(pi)
                equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
            
            equation_str = f"H = - [ " + " + ".join(equation_parts) + " ]"
            print(equation_str)

        print(f"Shannon's Diversity (H) = {shannon_index:.2f}\n")
        diversity_results.append(f"{shannon_index:.2f}")
    
    # Step 5: Phrase the final answer
    final_answer_string = ", ".join(diversity_results)
    print("Based on the calculations, the final answer is:")
    print(f"The insect diversity of each site is {final_answer_string}")
    print(f"<<<{final_answer_string}>>>")


if __name__ == '__main__':
    calculate_shannon_diversity()
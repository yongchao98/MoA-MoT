import math

def calculate_shannon_diversity(species_counts):
    """Calculates the Shannon's diversity index for a list of species counts."""
    total_individuals = sum(species_counts)
    
    if total_individuals == 0 or len(species_counts) <= 1:
        # Diversity is 0 if there's only one species or no individuals
        return 0.0, "H = -[ (1.0 * ln(1.0)) ]"

    shannon_index = 0
    equation_parts = []
    
    for count in species_counts:
        if count > 0:
            proportion = count / total_individuals
            shannon_index -= proportion * math.log(proportion)
            equation_parts.append(f"({proportion:.2f} * ln({proportion:.2f}))")
            
    equation_str = f"H = -[ {' + '.join(equation_parts)} ]"
    
    return shannon_index, equation_str

def solve_diversity_problem():
    """
    Solves the entomologist's diversity problem by identifying true species,
    recalculating abundances, and computing Shannon's diversity for each site.
    """
    # Step 1: Define the mapping from morphospecies to true species.
    # From visual inspection, A, B, C, D, and E are Harmonia axyridis. F is a spider.
    species_map = {
        'A': 'Harmonia axyridis',
        'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis',
        'E': 'Harmonia axyridis',
        'F': 'Spider'
    }

    # Step 2: Original abundance data based on morphospecies
    site_morpho_data = [
        {'A': 1, 'B': 1, 'C': 2},        # Site 1
        {'A': 10, 'D': 12, 'E': 1},     # Site 2
        {'B': 2, 'E': 2, 'F': 4}         # Site 3
    ]

    # Recalculate abundances based on true species for each site
    true_abundances_by_site = []
    for site_data in site_morpho_data:
        true_counts = {}
        for morpho, count in site_data.items():
            true_species = species_map[morpho]
            true_counts[true_species] = true_counts.get(true_species, 0) + count
        true_abundances_by_site.append(list(true_counts.values()))
    
    results = []
    print("Calculating Shannon's Diversity for each site based on true species:\n")

    # Step 3 & 4: Calculate and print diversity for each site
    for i, abundances in enumerate(true_abundances_by_site):
        site_num = i + 1
        diversity, equation = calculate_shannon_diversity(abundances)
        results.append(f"{diversity:.2f}")
        
        print(f"Site {site_num}:")
        print(f"  True species counts: {abundances}")
        print(f"  Calculation: {equation}")
        print(f"  Shannon's Diversity (H): {diversity:.2f}\n")
    
    # Step 5: Format and print the final answer
    final_answer_str = ", ".join(results)
    print("The final insect diversity of each site is:")
    print(final_answer_str)
    
    # Required final answer format
    print(f"<<<{final_answer_str}>>>")

# Run the solution
solve_diversity_problem()
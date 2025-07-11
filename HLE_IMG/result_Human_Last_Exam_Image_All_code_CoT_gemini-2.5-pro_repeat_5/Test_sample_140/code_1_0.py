import math

def calculate_shannon_diversity(species_counts, site_name):
    """Calculates the Shannon diversity index for a given list of species abundances."""
    total_individuals = sum(species_counts)
    
    # If there are no individuals or only one species, diversity is 0.
    if total_individuals == 0 or len(species_counts) <= 1:
        print(f"\nCalculating for {site_name}:")
        if total_individuals > 0:
            n1 = species_counts[0]
            p1 = n1 / total_individuals
            print(f"H' = - [ ({n1}/{total_individuals}) * ln({n1}/{total_individuals}) ]")
            print(f"H' = - [ {p1:.2f} * ln({p1:.2f}) ] = - [ {p1:.2f} * 0 ] = 0.00")
        else:
             print("H' = 0.00 (No individuals)")
        return 0.0

    shannon_index = 0
    equation_parts = []
    
    print(f"\nCalculating for {site_name}:")
    # Build the equation string first
    for count in species_counts:
        equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
    print(f"H' = - [ {' + '.join(equation_parts)} ]")
    
    # Now perform the calculation
    for count in species_counts:
        # p_i = proportion of individuals of species i
        p_i = count / total_individuals
        if p_i > 0:
            shannon_index -= p_i * math.log(p_i)
            
    print(f"H' = {shannon_index:.4f}")
    return shannon_index

def solve_diversity_problem():
    """
    Solves the entomology diversity problem based on true species identification.
    """
    # Step 1 & 2: Define data based on true species identification.
    # True Species 1: Harmonia axyridis (Morphospecies A, B, C, D, E)
    # True Species 2: Species F (Morphospecies F)
    
    # Abundances of true species at each site
    site1_abundances = [1 + 1 + 2]  # A + B + C = 4 of H. axyridis
    site2_abundances = [10 + 12 + 1] # A + D + E = 23 of H. axyridis
    site3_abundances = [2 + 2, 4]   # (B + E) = 4 of H. axyridis, F = 4 of Species F

    # Step 3: Calculate diversity for each site
    h1 = calculate_shannon_diversity(site1_abundances, "Site 1")
    h2 = calculate_shannon_diversity(site2_abundances, "Site 2")
    h3 = calculate_shannon_diversity(site3_abundances, "Site 3")
    
    # Step 4: Format and print the final result
    final_answer = f"{h1:.2f}, {h2:.2f}, {h3:.2f}"
    print("\nFinal Answer:")
    print(final_answer)

solve_diversity_problem()
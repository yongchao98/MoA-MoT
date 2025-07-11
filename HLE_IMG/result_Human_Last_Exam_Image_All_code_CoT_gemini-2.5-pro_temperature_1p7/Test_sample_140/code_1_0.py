import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon diversity for three insect survey sites after
    re-evaluating species identity based on true taxonomy.
    """
    
    # As identified from the images, morphospecies A-E are one species,
    # and F is a second, distinct species.

    # Aggregated counts based on true species identity
    sites_data = {
        'Site 1': {'Species 1': 1 + 1 + 2},
        'Site 2': {'Species 1': 10 + 12 + 1},
        'Site 3': {'Species 1': 2 + 2, 'Species 2': 4}
    }

    shannon_results = []
    
    print("Calculating Shannon's Diversity (H) for each site.")
    print("Formula: H = -Σ(pi * ln(pi)), where pi is the proportion of species i.\n")

    for site_name, site_counts in sites_data.items():
        print(f"--- {site_name} ---")
        
        total_individuals = sum(site_counts.values())
        shannon_index = 0.0
        
        equation_parts = []
        if total_individuals > 0:
            for species, count in site_counts.items():
                if count > 0:
                    p_i = count / total_individuals
                    # The Shannon index formula involves a sum over p_i * ln(p_i)
                    term = p_i * math.log(p_i)
                    shannon_index -= term
                    equation_parts.append(f"({count}/{total_individuals} * ln({count}/{total_individuals}))")
        
        # When a site has only one species, p_i is 1, and ln(1) is 0, making H = 0.
        if len(site_counts) <= 1:
            shannon_index = 0.0

        equation_str = "H = - [ " + " + ".join(equation_parts) + " ]"
        
        print(f"True Species Counts: {site_counts}")
        print(f"Total Individuals: {total_individuals}")
        print("Calculation:")
        print(f"{equation_str} = {shannon_index:.2f}")
        print("-" * (len(site_name) + 8) + "\n")
        
        shannon_results.append(shannon_index)
        
    final_answer = ", ".join([f"{h:.2f}" for h in shannon_results])
    print(f"The insect diversity (Shannon's H) for Site 1, Site 2, and Site 3 is:")
    print(final_answer)
    
    # Storing the final result for automated checking.
    global final_answer_value
    final_answer_value = final_answer

calculate_shannon_diversity()
<<<0.00, 0.00, 0.69>>>
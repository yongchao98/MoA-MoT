import math

def calculate_and_print_shannon_diversity():
    """
    Calculates and prints the Shannon's diversity for three insect survey sites
    after correcting for true species identity.
    """
    # Step 1: Define initial data and species mapping
    morpho_abundances = {
        'Site 1': {'A': 1, 'B': 1, 'C': 2},
        'Site 2': {'A': 10, 'D': 12, 'E': 1},
        'Site 3': {'B': 2, 'E': 2, 'F': 4}
    }

    species_map = {
        'A': 'Species 1 (H. axyridis)', 'B': 'Species 1 (H. axyridis)',
        'C': 'Species 1 (H. axyridis)', 'E': 'Species 1 (H. axyridis)',
        'D': 'Species 2',
        'F': 'Species 3'
    }

    print("Based on taxonomic identification, morphospecies A, B, C, and E are the same species.")
    print("The Shannon's diversity (H) is calculated as H = -Σ(pi * ln(pi)).\n")
    
    results = []
    # Step 2 & 3: Process each site
    for site_name in sorted(morpho_abundances.keys()):
        abundances = morpho_abundances[site_name]
        print(f"--- {site_name} ---")
        
        # Recalculate abundances based on true species
        true_species_counts = {}
        for morpho, count in abundances.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
        counts_list = list(true_species_counts.values())
        total_individuals = sum(counts_list)
        
        print(f"The true species abundances are {counts_list} for a total of {total_individuals} individuals.")

        shannon_val = 0.0
        
        # Handle cases with 0 or 1 species
        if total_individuals == 0 or len(counts_list) <= 1:
            p = 1.0
            print(f"The proportion (p1) is {counts_list[0]}/{total_individuals} = {p:.4f}")
            print(f"The equation is: H = -[({p:.4f} * ln({p:.4f}))]")
            print("Since ln(1) = 0, the diversity is 0.00.")
        else:
            equation_parts = []
            proportion_strings = []
            
            for i, count in enumerate(counts_list):
                proportion = count / total_individuals
                proportion_strings.append(f"p{i+1} = {count}/{total_individuals} = {proportion:.4f}")
                shannon_val -= proportion * math.log(proportion)
                # Format each number for the final equation output
                equation_parts.append(f"({proportion:.4f} * ln({proportion:.4f}))")
            
            print(f"The proportions are: {', '.join(proportion_strings)}")
            equation_str = "H = -[" + " + ".join(equation_parts) + "]"
            print(f"The equation is: {equation_str}")
            print(f"Shannon's Diversity (H) for {site_name} = {shannon_val:.2f}")

        results.append(f"{shannon_val:.2f}")
        print("-" * (len(site_name) + 6) + "\n")

    # Final combined answer
    final_answer = ", ".join(results)
    print(f"The insect diversity of each site to two decimal places is:\n{final_answer}")

calculate_and_print_shannon_diversity()
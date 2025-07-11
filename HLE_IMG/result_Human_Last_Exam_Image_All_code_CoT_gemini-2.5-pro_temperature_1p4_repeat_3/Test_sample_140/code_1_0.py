import math

def calculate_and_print_shannon_diversity():
    """
    Calculates and prints the Shannon diversity for three insect survey sites
    after correcting for true species identity.
    """
    
    # Original morphospecies abundance data
    sites_morpho = {
        1: {'A': 1, 'B': 1, 'C': 2},
        2: {'A': 10, 'D': 12, 'E': 1},
        3: {'B': 2, 'E': 2, 'F': 4}
    }
    
    # Mapping of morphospecies to true species
    species_map = {
        'A': 'H. axyridis', 
        'B': 'H. axyridis', 
        'C': 'H. axyridis',
        'D': 'H. axyridis', 
        'E': 'H. axyridis', 
        'F': 'Species_F'
    }
    
    # --- Data Consolidation ---
    sites_true_species = {}
    for site_num, morpho_counts in sites_morpho.items():
        true_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_counts[true_species] = true_counts.get(true_species, 0) + count
        sites_true_species[site_num] = true_counts
        
    final_diversities = []
    
    print("Calculating Shannon's Diversity for each site...\n")

    # --- Calculation for each site ---
    for site_num, true_counts in sites_true_species.items():
        print(f"--- Site {site_num} ---")
        print(f"True species counts: {true_counts}")
        
        abundances = list(true_counts.values())
        total_individuals = sum(abundances)
        
        print(f"Total individuals (N): {total_individuals}")
        
        shannon_index = 0.0
        
        if total_individuals > 0 and len(abundances) > 1:
            equation_parts = []
            proportions_info = []
            
            for species, count in true_counts.items():
                proportion = count / total_individuals
                proportions_info.append(f"p_{species} = {count} / {total_individuals} = {proportion:.2f}")
                term = proportion * math.log(proportion)
                shannon_index += term
                equation_parts.append(f"({proportion:.2f} * ln({proportion:.2f}))")
            
            # Print proportions and the full equation
            for p_info in proportions_info:
                print(p_info)
            print(f"Shannon's H = - ( {' + '.join(equation_parts)} )")
            
            shannon_index = -shannon_index

        elif len(abundances) == 1:
            # Diversity is 0 if only one species is present
            species_name = list(true_counts.keys())[0]
            count = list(true_counts.values())[0]
            proportion = count/total_individuals
            print(f"p_{species_name} = {count} / {total_individuals} = {proportion:.1f}")
            print(f"Shannon's H = - ({proportion:.1f} * ln({proportion:.1f}))")
            shannon_index = 0.0
        
        print(f"H = {shannon_index:.2f}\n")
        final_diversities.append(f"{shannon_index:.2f}")
        
    # --- Final Answer ---
    final_answer_string = ", ".join(final_diversities)
    print(f"The insect diversity of each site is: {final_answer_string}")
    print(f"<<<{final_answer_string}>>>")


calculate_and_print_shannon_diversity()
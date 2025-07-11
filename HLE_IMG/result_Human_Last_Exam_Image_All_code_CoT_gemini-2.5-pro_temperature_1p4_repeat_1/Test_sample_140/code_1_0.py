import math

def main():
    """
    Calculates and prints the Shannon's diversity for three insect survey sites
    after consolidating morphospecies into true species.
    """
    # Step 1: Define initial data
    # Abundance of each morphospecies at each site
    site_morpho_data = {
        "Site 1": {"A": 1, "B": 1, "C": 2},
        "Site 2": {"A": 10, "D": 12, "E": 1},
        "Site 3": {"B": 2, "E": 2, "F": 4},
    }

    # Identification of true species from morphospecies
    # Based on visual analysis, A, B, C, D, E are all morphs of Harmonia axyridis. 
    # F is a distinct species.
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis', 'F': 'Species F'
    }

    # Step 2: Consolidate data based on true species identity
    site_true_species_data = {}
    print("Consolidating morphospecies counts into true species counts:")
    for site, morpho_counts in site_morpho_data.items():
        consolidated_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            consolidated_counts[true_species] = consolidated_counts.get(true_species, 0) + count
        site_true_species_data[site] = consolidated_counts
        print(f"- {site}: {consolidated_counts}")
    print("-" * 30)

    # Step 3: Calculate and print Shannon's diversity for each site
    shannon_results = []
    site_names = sorted(site_true_species_data.keys())

    for site in site_names:
        counts_dict = site_true_species_data[site]
        counts = list(counts_dict.values())
        total_individuals = sum(counts)
        
        print(f"Calculation for {site}:")
        
        # Handle cases with 0 or 1 species
        if len(counts) <= 1:
            shannon_diversity = 0.0
            print(f"The site contains only one species. Shannon's diversity is 0.")
            print(f"H = {shannon_diversity:.2f}\n")
        else:
            shannon_diversity = 0.0
            equation_parts = []
            for count in counts:
                if count > 0:
                    proportion = count / total_individuals
                    shannon_diversity -= proportion * math.log(proportion)
                    equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
            
            equation_str = " + ".join(equation_parts)
            print(f"H = -[ {equation_str} ]")
            print(f"H = {shannon_diversity:.2f}\n")
            
        shannon_results.append(f"{shannon_diversity:.2f}")

    # Final formatted output
    final_answer = ", ".join(shannon_results)
    print("=" * 30)
    print(f"The insect diversity of each site is: {final_answer}")

if __name__ == "__main__":
    main()
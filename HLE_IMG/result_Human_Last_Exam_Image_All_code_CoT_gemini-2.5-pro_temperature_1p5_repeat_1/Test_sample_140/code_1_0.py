import math

def calculate_shannon_diversity_and_print():
    """
    This script calculates the Shannon's diversity for three insect survey sites
    after consolidating morphospecies into true species based on visual identification.
    """
    # Step 1: Define species mapping and site data
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
        'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis', 'F': 'Species F'
    }
    site_data = {
        'Site 1': {'A': 1, 'B': 1, 'C': 2},
        'Site 2': {'A': 10, 'D': 12, 'E': 1},
        'Site 3': {'B': 2, 'E': 2, 'F': 4}
    }

    # Step 2: Recalculate abundances for each site based on true species
    all_site_abundances = {}
    for site, morpho_counts in site_data.items():
        true_species_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        all_site_abundances[site] = true_species_counts

    # Step 3 & 4: Calculate diversity and print results
    results = []
    
    # Site 1
    s1_counts = all_site_abundances['Site 1']
    s1_n = list(s1_counts.values())[0]
    s1_N = s1_n
    s1_p = s1_n / s1_N
    h1 = - (s1_p * math.log(s1_p))
    results.append(h1)
    print("Site 1: 1 species with 4 individuals.")
    print(f"H = - [ ({s1_n}/{s1_N} * ln({s1_n}/{s1_N})) ] = {h1:.2f}")
    print("-" * 20)

    # Site 2
    s2_counts = all_site_abundances['Site 2']
    s2_n = list(s2_counts.values())[0]
    s2_N = s2_n
    s2_p = s2_n / s2_N
    h2 = - (s2_p * math.log(s2_p))
    results.append(h2)
    print("Site 2: 1 species with 23 individuals.")
    print(f"H = - [ ({s2_n}/{s2_N} * ln({s2_n}/{s2_N})) ] = {h2:.2f}")
    print("-" * 20)

    # Site 3
    s3_counts = all_site_abundances['Site 3']
    s3_n1 = s3_counts['Harmonia axyridis']
    s3_n2 = s3_counts['Species F']
    s3_N = s3_n1 + s3_n2
    s3_p1 = s3_n1 / s3_N
    s3_p2 = s3_n2 / s3_N
    h3 = - ((s3_p1 * math.log(s3_p1)) + (s3_p2 * math.log(s3_p2)))
    results.append(h3)
    print("Site 3: 2 species, with 4 individuals each.")
    print(f"H = - [ ({s3_n1}/{s3_N} * ln({s3_n1}/{s3_N})) + ({s3_n2}/{s3_N} * ln({s3_n2}/{s3_N})) ] = {h3:.2f}")
    print("-" * 20)
    
    # Final formatted answer
    final_answer_string = ", ".join([f"{r:.2f}" for r in results])
    print(f"Final Answer: {final_answer_string}")

calculate_shannon_diversity_and_print()
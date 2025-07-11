import math

def calculate_shannon_diversity():
    """
    Calculates the Shannon's diversity for three insect survey sites
    after correcting for true species identity.
    """
    
    # Step 1: Define the mapping from morphospecies to true species.
    # Based on visual identification, A, B, C, D, and E are all color morphs
    # of the same species (Harmonia axyridis), while F is a different species.
    species_map = {
        'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis',
        'C': 'Harmonia axyridis', 'D': 'Harmonia axyridis',
        'E': 'Harmonia axyridis', 'F': 'Other Species'
    }

    # Original morphospecies abundance data for each site.
    sites_morpho_data = {
        'Site 1': {'A': 1, 'B': 1, 'C': 2},
        'Site 2': {'A': 10, 'D': 12, 'E': 1},
        'Site 3': {'B': 2, 'E': 2, 'F': 4}
    }

    # Step 2: Recalculate abundances based on true species identity.
    sites_true_abundances = {}
    for site, morpho_counts in sites_morpho_data.items():
        true_counts = {}
        for morpho, count in morpho_counts.items():
            true_species = species_map[morpho]
            true_counts[true_species] = true_counts.get(true_species, 0) + count
        sites_true_abundances[site] = true_counts
    
    print("Shannon's Diversity Calculation (H = -Σ(pi * ln(pi))):\n")
    
    final_results = []

    # --- Site 1 Calculation ---
    site1_counts = list(sites_true_abundances['Site 1'].values())
    total_s1 = sum(site1_counts)
    p1 = site1_counts[0] / total_s1
    h1 = -(p1 * math.log(p1))
    final_results.append(h1)
    
    print("Site 1: Contains 1 true species with 4 individuals.")
    print(f"H = -[({site1_counts[0]}/{total_s1}) * ln({site1_counts[0]}/{total_s1})]")
    print(f"H = -[1.0 * ln(1.0)] = {h1:.2f}\n")

    # --- Site 2 Calculation ---
    site2_counts = list(sites_true_abundances['Site 2'].values())
    total_s2 = sum(site2_counts)
    p1 = site2_counts[0] / total_s2
    h2 = -(p1 * math.log(p1))
    final_results.append(h2)

    print("Site 2: Contains 1 true species with 23 individuals.")
    print(f"H = -[({site2_counts[0]}/{total_s2}) * ln({site2_counts[0]}/{total_s2})]")
    print(f"H = -[1.0 * ln(1.0)] = {h2:.2f}\n")

    # --- Site 3 Calculation ---
    site3_counts = list(sites_true_abundances['Site 3'].values())
    total_s3 = sum(site3_counts)
    p1 = site3_counts[0] / total_s3
    p2 = site3_counts[1] / total_s3
    h3 = -((p1 * math.log(p1)) + (p2 * math.log(p2)))
    final_results.append(h3)

    print("Site 3: Contains 2 true species with 4 individuals each.")
    print(f"H = -[ (({site3_counts[0]}/{total_s3}) * ln({site3_counts[0]}/{total_s3})) + (({site3_counts[1]}/{total_s3}) * ln({site3_counts[1]}/{total_s3})) ]")
    print(f"H = -[ ({p1} * ln({p1})) + ({p2} * ln({p2})) ] = {h3:.2f}\n")

    # Step 4: Print the final formatted answer.
    final_answer_string = ", ".join([f"{res:.2f}" for res in final_results])
    print("The insect diversity of each site is:")
    print(final_answer_string)


calculate_shannon_diversity()
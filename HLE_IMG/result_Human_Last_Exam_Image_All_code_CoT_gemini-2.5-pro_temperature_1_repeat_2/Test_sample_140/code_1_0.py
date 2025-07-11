import math

def calculate_shannon(site_name, counts):
    """
    Calculates and prints the Shannon's diversity index for a given site.
    
    Args:
        site_name (str): The name of the site.
        counts (dict): A dictionary with species as keys and abundances as values.
    """
    N = sum(counts.values())
    
    print(f"\n--- {site_name} ---")
    print(f"Species Abundances: {counts}")
    print(f"Total Individuals (N): {N}")

    if N == 0 or len(counts) <= 1:
        # Diversity is 0 if there is only one species or no individuals
        shannon_index = 0.0
        species_name = list(counts.keys())[0] if counts else 'N/A'
        n_i = N
        p_i = 1.0 if N > 0 else 0.0
        print(f"Calculation: H' = - (({n_i}/{N}) * ln({n_i}/{N}))")
    else:
        shannon_index = 0.0
        equation_parts = []
        for species, n_i in counts.items():
            if n_i > 0:
                p_i = n_i / N
                shannon_index -= p_i * math.log(p_i)
                equation_parts.append(f"({n_i}/{N}) * ln({n_i}/{N})")
        
        equation_str = " + ".join(equation_parts)
        print(f"Calculation: H' = - ({equation_str})")

    print(f"Shannon's Diversity (H'): {shannon_index:.2f}")
    return shannon_index

# 1. Define the initial survey data and the taxonomic mapping
site_data = {
    "Site 1": {"A": 1, "B": 1, "C": 2},
    "Site 2": {"A": 10, "D": 12, "E": 1},
    "Site 3": {"B": 2, "E": 2, "F": 4}
}

# A, B, C, D, E are the same species (S1: Harmonia axyridis)
# F is a different species (S2)
morpho_to_true_species = {
    'A': 'S1', 'B': 'S1', 'C': 'S1',
    'D': 'S1', 'E': 'S1', 'F': 'S2'
}

print("Based on taxonomic analysis, morphospecies A, B, C, D, and E are identified as color variants of the same species (S1), while F is a distinct species (S2).")
print("Recalculating diversity for each site based on true species counts:")

site_diversities = []

# 2. Loop through each site to recalculate abundances and diversity
for site_name, morpho_counts in site_data.items():
    true_species_counts = {}
    for morpho, count in morpho_counts.items():
        true_species = morpho_to_true_species[morpho]
        # Add the count to the correct true species
        true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
    # 3. Calculate Shannon's diversity for the site
    diversity = calculate_shannon(site_name, true_species_counts)
    site_diversities.append(diversity)

# 4. Format and print the final answer
final_answer = ", ".join([f"{d:.2f}" for d in site_diversities])
print("\nFinal Answer:")
print(final_answer)

<<<0.00, 0.00, 0.69>>>
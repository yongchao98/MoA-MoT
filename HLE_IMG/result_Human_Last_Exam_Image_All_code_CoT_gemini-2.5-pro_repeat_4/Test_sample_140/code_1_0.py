import math

def calculate_and_print_shannon_diversity(site_name, abundances):
    """
    Calculates and prints the Shannon diversity index and its components.
    
    Args:
        site_name (str): The name of the site for printing.
        abundances (dict): A dictionary with species names as keys and counts as values.
    
    Returns:
        float: The calculated Shannon diversity index.
    """
    print(f"--- {site_name} ---")
    
    counts = list(abundances.values())
    total_individuals = sum(counts)
    
    print(f"True Species Abundances: {abundances}")
    print(f"Total individuals (N): {total_individuals}")
    
    # If there's 0 or 1 species, diversity is 0.
    if total_individuals == 0 or len(counts) <= 1:
        print("Shannon's Diversity (H) = 0.00 (as there is only one species type).")
        print("-" * 25)
        return 0.0
        
    shannon_index = 0.0
    equation_parts = []
    
    # Build the equation string and calculate the index
    for species, count in abundances.items():
        if count > 0:
            proportion = count / total_individuals
            shannon_index -= proportion * math.log(proportion)
            equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
    
    print("Shannon's Diversity H = - [ " + " + ".join(equation_parts) + " ]")
    print(f"H = {shannon_index:.2f}")
    print("-" * 25)
    
    return shannon_index

# --- Step 1: Define initial data and species mapping ---

# Morphospecies abundances from the survey
site1_morpho = {'A': 1, 'B': 1, 'C': 2}
site2_morpho = {'A': 10, 'D': 12, 'E': 1}
site3_morpho = {'B': 2, 'E': 2, 'F': 4}

# Mapping from morphospecies to true species
# A, B, C, D, E are all Harmonia axyridis. F is a different species.
morpho_to_species = {
    'A': 'Harmonia axyridis',
    'B': 'Harmonia axyridis',
    'C': 'Harmonia axyridis',
    'D': 'Harmonia axyridis',
    'E': 'Harmonia axyridis',
    'F': 'Species F' 
}

# --- Step 2: Update site abundances with true species counts ---

all_sites_data = [site1_morpho, site2_morpho, site3_morpho]
true_species_by_site = []

for site_data in all_sites_data:
    aggregated_counts = {}
    for morpho, count in site_data.items():
        true_species = morpho_to_species[morpho]
        aggregated_counts[true_species] = aggregated_counts.get(true_species, 0) + count
    true_species_by_site.append(aggregated_counts)
    
# --- Step 3: Calculate Shannon's diversity for each site ---

results = []
for i, site_abundances in enumerate(true_species_by_site):
    site_name = f"Site {i + 1}"
    diversity = calculate_and_print_shannon_diversity(site_name, site_abundances)
    results.append(diversity)

# --- Step 4: Format and print the final answer ---

final_answer_string = ", ".join([f"{r:.2f}" for r in results])
print("\nFinal Answer:")
print(f"The insect diversity of each site is: {final_answer_string}")
<<<0.00, 0.00, 0.69>>>
import math

def calculate_and_print_shannon(site_name, true_species_counts):
    """
    Calculates and prints the Shannon diversity index and the equation used.
    
    Args:
        site_name (str): The name of the site (e.g., "Site 1").
        true_species_counts (dict): A dictionary with species names as keys and counts as values.
    """
    print(f"--- {site_name} ---")
    counts = list(true_species_counts.values())
    total_individuals = sum(counts)
    
    print(f"True Species Counts: {true_species_counts}")
    
    if total_individuals == 0:
        print("No individuals found, diversity is 0.")
        return 0.0

    shannon_h = 0.0
    equation_parts = []
    
    # Build the equation string and calculate H
    for count in counts:
        proportion = count / total_individuals
        if proportion > 0:
            shannon_h -= proportion * math.log(proportion)
            equation_parts.append(f"({proportion:.4f} * ln({proportion:.4f}))")

    # Final equation format
    # The negative sign is applied to the sum of all parts
    equation_string = f"H = -[ {' + '.join(equation_parts)} ]"
    
    print("Equation:")
    print(f"{equation_string} = {shannon_h:.2f}\n")
    return shannon_h

# Original abundance data by morphospecies
site1_morpho = {'A': 1, 'B': 1, 'C': 2}
site2_morpho = {'A': 10, 'D': 12, 'E': 1}
site3_morpho = {'B': 2, 'E': 2, 'F': 4}

# Mapping of morphospecies to true species
# True Species 1: Harmonia axyridis (A, B, C, D, E)
# True Species 2: Distinct species (F)
true_species_map = {
    'A': 'Harmonia axyridis', 'B': 'Harmonia axyridis', 'C': 'Harmonia axyridis',
    'D': 'Harmonia axyridis', 'E': 'Harmonia axyridis', 'F': 'Species F'
}

all_sites_morpho = {
    "Site 1": site1_morpho,
    "Site 2": site2_morpho,
    "Site 3": site3_morpho
}

# Recalculate abundances based on true species and compute diversity
results = []
for site_name, morpho_counts in all_sites_morpho.items():
    true_counts = {}
    for morpho, count in morpho_counts.items():
        true_spec = true_species_map[morpho]
        true_counts[true_spec] = true_counts.get(true_spec, 0) + count
    
    h_value = calculate_and_print_shannon(site_name, true_counts)
    results.append(h_value)

# Format and print the final answer
final_answer_string = ", ".join([f"{h:.2f}" for h in results])
print("The insect diversity of each site is:")
print(final_answer_string)
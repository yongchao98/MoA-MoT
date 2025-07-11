import math

def calculate_shannon_and_print(site_name, true_species_counts):
    """
    Calculates Shannon's diversity and prints a step-by-step explanation.
    """
    print(f"----- {site_name} Analysis -----")
    print(f"True species counts: {true_species_counts}")

    total_individuals = sum(true_species_counts.values())
    if total_individuals == 0:
        print("No individuals found. Diversity is 0.00.")
        return 0.0
        
    print(f"Total individuals (N): {total_individuals}")
    
    # Check for sites with only one species
    if len(true_species_counts) <= 1:
        p = 1.0
        print(f"Only one species present. The proportion (p1) is {list(true_species_counts.values())[0]} / {total_individuals} = {p:.1f}.")
        print(f"Equation: H' = -[({p:.1f}) * ln({p:.1f})]")
        print("Since ln(1) = 0, the Shannon's Diversity (H') is 0.00.\n")
        return 0.0

    proportions = [count / total_individuals for count in true_species_counts.values()]
    
    # Build and print the equation parts
    equation_parts = []
    for p in proportions:
        equation_parts.append(f"({p:.2f} * ln({p:.2f}))")
    
    print("Shannon's Diversity Equation: H' = -[ " + " + ".join(equation_parts) + " ]")

    # Calculate diversity
    shannon_index = -sum(p * math.log(p) for p in proportions)
    
    print(f"Calculated Shannon's Diversity (H'): {shannon_index:.2f}\n")
    return shannon_index

# Step 1: Define morphospecies to true species mapping
# A, B, C, D, E -> Harmonia axyridis (True Species 1)
# F -> Different species (True Species 2)
morpho_to_true = {
    'A': 'True Species 1', 'B': 'True Species 1', 'C': 'True Species 1',
    'D': 'True Species 1', 'E': 'True Species 1', 'F': 'True Species 2'
}

# Step 2: Define morphospecies abundance data for each site
sites_morpho_data = {
    "Site 1": {'A': 1, 'B': 1, 'C': 2},
    "Site 2": {'A': 10, 'D': 12, 'E': 1},
    "Site 3": {'B': 2, 'E': 2, 'F': 4}
}

# Step 3: Recalculate abundances and calculate diversity for each site
site_diversities = []
for site_name, morpho_counts in sites_morpho_data.items():
    true_counts = {}
    for morpho_species, count in morpho_counts.items():
        true_species_id = morpho_to_true[morpho_species]
        true_counts[true_species_id] = true_counts.get(true_species_id, 0) + count
    
    diversity = calculate_shannon_and_print(site_name, true_counts)
    site_diversities.append(diversity)

# Final formatted output
final_answer = ", ".join([f"{d:.2f}" for d in site_diversities])
print("--------------------------")
print(f"The final calculated insect diversity of each site is: {final_answer}")
print("<<<" + final_answer + ">>>")
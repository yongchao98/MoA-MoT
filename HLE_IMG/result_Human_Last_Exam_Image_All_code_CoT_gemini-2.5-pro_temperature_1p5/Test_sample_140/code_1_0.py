import math

def calculate_shannon_with_explanation(counts_dict, site_name):
    """
    Calculates Shannon's diversity index for a given dictionary of species counts
    and prints the formula with the specific numbers used.
    """
    
    # Get the counts of individuals for each species
    counts = list(counts_dict.values())
    # Filter out any species with 0 individuals
    counts = [c for c in counts if c > 0]
    
    total_individuals = sum(counts)
    
    print(f"Calculating Shannon's diversity for {site_name}:")

    # Handle cases with 0 or 1 species, where diversity is 0.
    if total_individuals == 0:
        print("No insects found. Total individuals (N) = 0.")
        print("H' = 0.00")
        return 0.0
    if len(counts) <= 1:
        n1 = counts[0]
        p1 = n1 / total_individuals
        # The equation for a single species
        print(f"Only one species present. Total individuals (N) = {total_individuals}.")
        print(f"H' = -[({n1}/{total_individuals}) * ln({n1}/{total_individuals})]")
        # ln(1) is 0, so the result is 0
        print("= 0.00")
        return 0.0

    shannon_sum = 0.0
    equation_parts = []
    
    # Build the equation string and calculate the sum for multiple species
    for count in counts:
        proportion = count / total_individuals
        shannon_sum += proportion * math.log(proportion)
        equation_parts.append(f"({count}/{total_individuals} * ln({count}/{total_individuals}))")

    equation_str = f"H' = -[ {' + '.join(equation_parts)} ]"
    final_diversity = -shannon_sum
    
    print(f"Total individuals (N) = {total_individuals}.")
    print(equation_str)
    print(f"= {final_diversity:.2f}")
    
    return final_diversity

# Step 1: Map morphospecies to their true taxonomic group.
# A, B, C, E are morphs of Harmonia axyridis.
# D is a distinct lady beetle species.
# F is a spider (not an insect) and will be excluded.
morpho_to_true_species = {
    'A': 'Harmonia axyridis',
    'B': 'Harmonia axyridis',
    'C': 'Harmonia axyridis',
    'D': 'Species D',
    'E': 'Harmonia axyridis',
    'F': 'Spider (Not an Insect)'
}

# Original abundance data
site_abundances = [
    {'A': 1, 'B': 1, 'C': 2},        # Site 1
    {'A': 10, 'D': 12, 'E': 1},      # Site 2
    {'B': 2, 'E': 2, 'F': 4}         # Site 3
]

shannon_results = []
# Step 2 & 3: Consolidate data and calculate diversity for each site
for i, site_data in enumerate(site_abundances):
    true_species_counts = {}
    for morpho, count in site_data.items():
        true_species = morpho_to_true_species[morpho]
        
        # Exclude non-insects from the calculation
        if "Not an Insect" in true_species:
            continue
            
        # Sum the counts for the true species
        true_species_counts[true_species] = true_species_counts.get(true_species, 0) + count
        
    diversity = calculate_shannon_with_explanation(true_species_counts, f"Site {i+1}")
    shannon_results.append(diversity)
    print("-" * 30)

# Step 4: Format and print the final answer
final_answer_str = ", ".join([f"{result:.2f}" for result in shannon_results])
print("The insect diversity of each site is:")
print(final_answer_str)
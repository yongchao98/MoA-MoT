import math

def calculate_shannon_diversity(species_counts):
    """
    Calculates the Shannon's diversity index (H').
    The formula is H' = -Σ(pi * ln(pi)), where pi is the proportion of species i.
    """
    total_individuals = sum(species_counts)
    if total_individuals == 0:
        return 0.0

    shannon_index = 0.0
    for count in species_counts:
        if count > 0:
            proportion = count / total_individuals
            shannon_index += proportion * math.log(proportion)

    return -shannon_index

# Step 1: Define the raw abundance data from the survey
site1_morpho = {'A': 1, 'B': 1, 'C': 2}
site2_morpho = {'A': 10, 'D': 12, 'E': 1}
site3_morpho = {'B': 2, 'E': 2, 'F': 4}

# Step 2: Define the true species based on the taxonomist's identification
# Morphospecies A, B, C, and E are the same species (Harmonia axyridis).
# D and F are distinct species.
species_map = {
    'A': 'Harmonia axyridis',
    'B': 'Harmonia axyridis',
    'C': 'Harmonia axyridis',
    'D': 'Species D',
    'E': 'Harmonia axyridis',
    'F': 'Species F'
}

sites_data = [site1_morpho, site2_morpho, site3_morpho]
final_shannon_values = []

print("Calculating Shannon's Diversity for each site after taxonomic correction:\n")

# Step 3 & 4: Recalculate abundances and compute diversity for each site
for i, site_data in enumerate(sites_data):
    site_num = i + 1
    true_abundances = {}
    # Recalculate abundances based on true species
    for morpho, count in site_data.items():
        true_species = species_map[morpho]
        true_abundances[true_species] = true_abundances.get(true_species, 0) + count

    counts = list(true_abundances.values())
    total = sum(counts)
    
    print(f"--- Site {site_num} ---")
    print(f"Original Morphospecies Counts: {site_data}")
    print(f"Corrected True Species Counts: {true_abundances}")

    # Build and print the equation string
    equation_parts = []
    for count in counts:
        p = count / total
        equation_parts.append(f"({p:.3f} * ln({p:.3f}))")
    
    equation_str = "H' = -[" + " + ".join(equation_parts) + "]"
    print(f"Calculation: {equation_str}")
    
    shannon_value = calculate_shannon_diversity(counts)
    final_shannon_values.append(shannon_value)
    print(f"Shannon's Diversity (H') = {shannon_value:.2f}\n")

# Final formatted output
final_answer_str = ", ".join([f"{val:.2f}" for val in final_shannon_values])
print("--- Final Answer ---")
print(f"The Shannon's diversity of insects at each site is: {final_answer_str}")

# The final answer in the required format
print(f"\n<<< {final_answer_str} >>>")
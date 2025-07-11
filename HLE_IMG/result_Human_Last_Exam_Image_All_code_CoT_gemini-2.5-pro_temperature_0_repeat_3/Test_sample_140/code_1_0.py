import math

def calculate_shannon_diversity(site_name, abundances):
    """
    Calculates and prints the Shannon's diversity index for a given site.
    
    Args:
        site_name (str): The name of the site (e.g., "Site 1").
        abundances (list): A list of abundances for each species at the site.
    """
    total_individuals = sum(abundances)
    
    if total_individuals == 0:
        print(f"{site_name}: H = 0.00")
        return 0.00

    shannon_index = 0
    equation_parts = []
    
    for count in abundances:
        if count > 0:
            proportion = count / total_individuals
            shannon_index -= proportion * math.log(proportion)
            equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")

    equation_str = " + ".join(equation_parts)
    print(f"{site_name}: H = - ( {equation_str} ) = {shannon_index:.2f}")
    return round(shannon_index, 2)

# Step 1 & 2: Consolidate abundances based on true species identification.
# True Species 1: A, B, C, E
# True Species 2: D
# True Species 3: F

# Original Data:
# Site 1) A: 1, B: 1, C: 2
# Site 2) A: 10, D: 12, E: 1
# Site 3) B: 2, E: 2, F: 4

# New abundances based on true species:
site1_abundances = [1 + 1 + 2]  # All are True Species 1
site2_abundances = [10 + 1, 12] # True Species 1, True Species 2
site3_abundances = [2 + 2, 4]   # True Species 1, True Species 3

# Step 3: Calculate Shannon's diversity for each site and print the steps.
print("Calculating Shannon's Diversity (H) for each site:")
h1 = calculate_shannon_diversity("Site 1", site1_abundances)
h2 = calculate_shannon_diversity("Site 2", site2_abundances)
h3 = calculate_shannon_diversity("Site 3", site3_abundances)

# Step 4: Format and print the final answer.
final_answer = f"{h1:.2f}, {h2:.2f}, {h3:.2f}"
print("\nFinal Answer:")
print(final_answer)
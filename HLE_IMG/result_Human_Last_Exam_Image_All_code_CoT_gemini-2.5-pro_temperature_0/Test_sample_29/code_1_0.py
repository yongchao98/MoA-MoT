import numpy as np

def calculate_and_print_hill_number(level_name, abundances, q):
    """
    Calculates the Hill number for a given taxonomic level and prints the steps.
    
    Args:
        level_name (str): The name of the taxonomic level (e.g., "Species").
        abundances (list): A list of counts for each group in the level.
        q (int): The order of diversity 'q'.
        
    Returns:
        float: The calculated Hill number.
    """
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    
    # The numbers (abundances) used in the equation
    print(f"Abundances: {abundances}")
    
    total_individuals = sum(abundances)
    if total_individuals == 0:
        print("No individuals found. Hill Number is 0.\n")
        return 0.0
        
    proportions = np.array(abundances) / total_individuals
    # Filter out zero proportions to avoid issues with log(0)
    proportions = proportions[proportions > 0]
    
    # The proportions (p_i) used in the equation
    print(f"Proportions (p_i): {[f'{p:.4f}' for p in proportions]}")

    if q == 1:
        # For q=1, Hill number is the exponential of Shannon entropy
        shannon_entropy = -np.sum(proportions * np.log(proportions))
        print(f"Shannon Entropy (-Σ p_i * ln(p_i)): {shannon_entropy:.4f}")
        hill_val = np.exp(shannon_entropy)
        print(f"Final Equation: D_1 = exp({shannon_entropy:.4f})")
    else:
        # For q != 1, use the general formula
        sum_p_q = np.sum(proportions ** q)
        print(f"Sum of powered proportions (Σ p_i^{q}): {sum_p_q:.4f}")
        exponent = 1 / (1 - q)
        hill_val = sum_p_q ** exponent
        print(f"Final Equation: D_{q} = ({sum_p_q:.4f}) ^ (1 / (1 - {q}))")

    print(f"Resulting Hill Number for {level_name}: {hill_val:.4f}\n")
    return hill_val

# Step 1 & 2: Define abundances for each taxonomic level based on visual inspection
# Species counts: 4 distinct species identified
abundances_species = [76, 28, 6, 8]

# Genus counts: Species 1 & 2 are grouped into one genus, Species 3 & 4 are in separate genera
abundances_genus = [76 + 28, 6, 8]

# Family counts: All specimens are assumed to be in one family
total_individuals = sum(abundances_species)
abundances_family = [total_individuals]

# Order counts: All specimens are in one order
abundances_order = [total_individuals]

# Step 3 & 4: Calculate Hill numbers for each level with its respective q
# q=1 for Order
hill_order = calculate_and_print_hill_number("Order", abundances_order, 1)

# q=2 for Family
hill_family = calculate_and_print_hill_number("Family", abundances_family, 2)

# q=3 for Genus
hill_genus = calculate_and_print_hill_number("Genus", abundances_genus, 3)

# q=4 for Species
hill_species = calculate_and_print_hill_number("Species", abundances_species, 4)

# Final formatted output
final_results = [hill_order, hill_family, hill_genus, hill_species]
formatted_results = ", ".join([f"{res:.2f}" for res in final_results])
print("--- Final Answer ---")
print("The Hill numbers for order (q=1), family (q=2), genus (q=3), and species (q=4) are:")
print(formatted_results)
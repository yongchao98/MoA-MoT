import numpy as np

def calculate_and_print_hill_number(level_name, abundances, q):
    """
    Calculates the Hill number and prints the equation components and result.
    
    Args:
        level_name (str): The name of the taxonomic level.
        abundances (list): A list of counts for each group.
        q (int): The order of diversity 'q'.

    Returns:
        float: The calculated Hill number.
    """
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    
    abundances = np.array(abundances, dtype=float)
    abundances = abundances[abundances > 0] # Exclude groups with zero individuals
    
    total_individuals = np.sum(abundances)
    proportions = abundances / total_individuals

    print(f"Abundances: {abundances.astype(int).tolist()}")
    print(f"Total Individuals: {int(total_individuals)}")
    
    # Handle the special case where q = 1 (Shannon diversity)
    if q == 1:
        # H' = -sum(p_i * ln(p_i))
        # D_1 = exp(H')
        shannon_entropy_components = proportions * np.log(proportions)
        shannon_entropy = -np.sum(shannon_entropy_components)
        result = np.exp(shannon_entropy)
        
        # Format the equation string
        equation_parts = [f"{p:.2f}*ln({p:.2f})" for p in proportions]
        equation = f"D_{q} = exp( - ( {' + '.join(equation_parts)} ) )"
        print(equation)

    # Handle the general case for q != 1
    else:
        # D_q = (sum(p_i^q))^(1 / (1 - q))
        sum_p_q = np.sum(proportions**q)
        exponent = 1 / (1 - q)
        result = sum_p_q**exponent

        # Format the equation string
        equation_parts = [f"({int(n)}/{int(total_individuals)})^{q}" for n in abundances]
        equation = f"D_{q} = ( {' + '.join(equation_parts)} ) ^ (1 / (1 - {q}))"
        print(equation)

    print(f"Result: {result:.4f}\n")
    return result

# --- Data from Image Analysis ---
# Taxonomic Level: Order (Lepidoptera)
# All 88 specimens are in one order.
order_abundances = [88]

# Taxonomic Level: Family (Nymphalidae)
# All 88 specimens are in one family.
family_abundances = [88]

# Taxonomic Level: Genus (Adelpha)
# All 88 specimens are in one genus.
genus_abundances = [88]

# Taxonomic Level: Species
# Three distinct species were identified with the following counts.
species_abundances = [48, 20, 20]

# --- Calculations ---
hill_order = calculate_and_print_hill_number("Order", order_abundances, 1)
hill_family = calculate_and_print_hill_number("Family", family_abundances, 2)
hill_genus = calculate_and_print_hill_number("Genus", genus_abundances, 3)
hill_species = calculate_and_print_hill_number("Species", species_abundances, 4)

# --- Final Output ---
final_results = [hill_order, hill_family, hill_genus, hill_species]
formatted_results = ", ".join([f"{res:.2f}" for res in final_results])
print("--- Final Answer ---")
print("The Hill numbers for q=1 (Order), q=2 (Family), q=3 (Genus), and q=4 (Species) are:")
print(formatted_results)

# Return the answer in the specified format
print(f"\n<<<{formatted_results}>>>")
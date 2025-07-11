import numpy as np

def calculate_and_print_hill(level_name, abundances, q):
    """
    Calculates and prints the Hill number for a given set of abundances and order q.
    """
    N = sum(abundances)
    # Filter out zero-abundance groups to avoid division by zero
    proportions = [n / N for n in abundances if n > 0]
    
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    print(f"Abundances (n_i): {abundances}, Total (N): {N}")
    
    # Handle the special case for q=1
    if q == 1:
        # D_1 = exp(H') where H' = -Σ(p_i * ln(p_i))
        shannon_entropy_terms = []
        shannon_equation_parts = []
        for p in proportions:
            shannon_entropy_terms.append(p * np.log(p))
            shannon_equation_parts.append(f"({p:.2f}*ln({p:.2f}))")
        
        shannon_entropy = -sum(shannon_entropy_terms)
        hill_value = np.exp(shannon_entropy)
        
        print(f"Formula: D_1 = exp( - ( { ' + '.join(shannon_equation_parts)} ) )")
        print(f"Resulting Hill Number ({level_name}): {hill_value:.2f}\n")

    # Handle the general case for q != 1
    else:
        # D_q = (Σ(p_i^q))^(1/(1-q))
        sum_of_powers = sum([p**q for p in proportions])
        exponent = 1 / (1 - q)
        hill_value = sum_of_powers**exponent
        
        power_equation_parts = [f"{p:.2f}^{q}" for p in proportions]
        print(f"Formula: D_{q} = ( { ' + '.join(power_equation_parts)} ) ^ (1/(1-{q}))")
        print(f"Resulting Hill Number ({level_name}): {hill_value:.2f}\n")
        
    return hill_value

# Based on visual analysis of the provided image:
# Total individuals N = 12+12+8+6+8+12+12+8+5+10 = 93
order_abundances = [93]         # All are Order Lepidoptera
family_abundances = [93]        # All appear to be Family Nymphalidae
genus_abundances = [75, 18]       # Grouped into 2 genera based on major pattern differences
species_abundances = [48, 16, 6, 18, 5] # Grouped into 5 species based on fine pattern differences

# Calculate the required Hill numbers
# q=1 for order, q=2 for family, q=3 for genus, q=4 for species
hill_order = calculate_and_print_hill("Order", order_abundances, 1)
hill_family = calculate_and_print_hill("Family", family_abundances, 2)
hill_genus = calculate_and_print_hill("Genus", genus_abundances, 3)
hill_species = calculate_and_print_hill("Species", species_abundances, 4)

# Print the final result in the requested format
final_answer = f"{hill_order:.2f},{hill_family:.2f},{hill_genus:.2f},{hill_species:.2f}"
print("---")
print("Final comma-separated values for Order (q=1), Family (q=2), Genus (q=3), Species (q=4):")
print(final_answer)
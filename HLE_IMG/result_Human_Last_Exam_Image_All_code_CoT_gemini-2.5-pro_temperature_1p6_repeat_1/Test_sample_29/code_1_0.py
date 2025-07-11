import numpy as np

def calculate_and_print_hill(counts, q, level_name):
    """
    Calculates the Hill number and prints the detailed calculation process.
    """
    counts = np.array(counts)
    N = np.sum(counts)
    # Filter out zero counts to avoid errors in log
    non_zero_counts = counts[counts > 0]
    pi = non_zero_counts / N

    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    print(f"Abundance Counts: {counts.tolist()}")
    print(f"Total Individuals (N): {N}")
    
    if q == 1:
        # Special case for q=1: exponential of Shannon entropy
        shannon_entropy = -np.sum(pi * np.log(pi))
        hill_number = np.exp(shannon_entropy)
        
        # Build and print the equation string
        equation_parts = [f"({p:.3f} * ln({p:.3f}))" for p in pi]
        print(f"Equation: D1 = exp(-[ {' + '.join(equation_parts)} ])")

    else:
        # General case for q != 1
        sum_pi_q = np.sum(pi**q)
        exponent = 1 / (1 - q)
        hill_number = sum_pi_q**exponent
        
        # Build and print the equation string
        equation_parts = [f"{p:.3f}^{q}" for p in pi]
        print(f"Equation: D{q} = ( {' + '.join(equation_parts)} )^(1/(1-{q}))")

    print(f"Result for {level_name}: {hill_number:.2f}")
    print("-" * 50)
    return hill_number

# Step 1: Counts derived from image analysis.
# Total individuals = 77 + 8 + 19 + 6 = 110
counts_species = [77, 8, 19, 6]

# Step 2: Group counts for higher taxonomic levels based on visual similarity.
# Genus/Family Level: Grouping visually similar species (77+8), (19), (6)
counts_genus_family = [85, 19, 6]
# Order Level: All specimens are in one order
counts_order = [110]

# Step 3: Calculate Hill numbers for the specified q values and taxonomic levels.
hill_order_q1 = calculate_and_print_hill(counts_order, 1, "Order")
hill_family_q2 = calculate_and_print_hill(counts_genus_family, 2, "Family")
hill_genus_q3 = calculate_and_print_hill(counts_genus_family, 3, "Genus")
hill_species_q4 = calculate_and_print_hill(counts_species, 4, "Species")

# Step 4: Print the final comma-separated result.
final_results = [hill_order_q1, hill_family_q2, hill_genus_q3, hill_species_q4]
formatted_results = ", ".join([f"{res:.2f}" for res in final_results])
print("\nFinal Result (Order, Family, Genus, Species):")
print(formatted_results)
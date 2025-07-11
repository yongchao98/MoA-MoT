import numpy as np

def calculate_hill(abundances, q, level_name):
    """Calculates and prints the Hill number for a given list of abundances and order q."""
    abundances_arr = np.array(abundances, dtype=float)
    total_individuals = np.sum(abundances_arr)
    
    if total_individuals == 0:
        print(f"Calculation for {level_name} (q={q}):")
        print("  - No individuals found. Diversity is 0.\n")
        return 0.0
        
    proportions = abundances_arr / total_individuals
    proportions_gt_zero = proportions[proportions > 0]
    
    print(f"Calculation for {level_name} (q={q}):")
    print(f"  - Abundances: {abundances}")
    print(f"  - Total Individuals: {int(total_individuals)}")
    print(f"  - Proportions (pi): {[round(p, 2) for p in proportions]}")
    
    # Special case for q = 1
    if q == 1:
        shannon_entropy = -np.sum(proportions_gt_zero * np.log(proportions_gt_zero))
        result = np.exp(shannon_entropy)
        print(f"  - Equation (q=1): D^1 = exp(-Σ(pi * ln(pi)))")
        print(f"  - D^1 = exp({shannon_entropy:.4f})")
    # General case for q != 1
    else:
        sum_p_q = np.sum(proportions_gt_zero**q)
        result = sum_p_q**(1 / (1 - q))
        print(f"  - Equation (q={q}): D^{q} = (Σ(pi^{q}))^(1/(1-{q}))")
        power_str = " + ".join([f"{p:.2f}^{q}" for p in proportions])
        print(f"  - D^{q} = ({power_str})^(1/{1-q})")
        print(f"  - D^{q} = ({sum_p_q:.4f})^({1/(1-q):.4f})")

    print(f"  - Result: {result:.4f}\n")
    return result

# --- Data based on analysis of the image ---

# Level 1: Order (q=1)
order_abundances = [100]
hill_order = calculate_hill(order_abundances, 1, "Order")

# Level 2: Family (q=2)
family_abundances = [100]
hill_family = calculate_hill(family_abundances, 2, "Family")

# Level 3: Genus (q=3)
genus_abundances = [90, 10]
hill_genus = calculate_hill(genus_abundances, 3, "Genus")

# Level 4: Species (q=4)
species_abundances = [58, 17, 15, 10]
hill_species = calculate_hill(species_abundances, 4, "Species")

# --- Final Answer ---
final_answer = f"{hill_order:.2f},{hill_family:.2f},{hill_genus:.2f},{hill_species:.2f}"
print("Final rounded values, separated by commas:")
print(final_answer)
print(f"<<<{final_answer}>>>")

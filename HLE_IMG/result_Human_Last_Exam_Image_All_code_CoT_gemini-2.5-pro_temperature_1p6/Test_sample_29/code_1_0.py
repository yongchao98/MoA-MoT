import numpy as np
import math

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates the Hill number for a given list of abundances and q value.
    It also prints the details of the calculation, including the formula.
    """
    # Convert to a numpy array for vectorized operations
    n = np.array(abundances, dtype=float)
    N = np.sum(n)
    
    # Avoid division by zero if there are no individuals
    if N == 0:
        p = np.zeros_like(n)
    else:
        p = n / N

    print(f"--- Calculation for {level_name} Level (q={q}) ---")
    print(f"Abundance counts (n_i): {abundances}")
    print(f"Total individuals (N): {int(N)}")
    
    # Generate strings for the formula parts
    proportions_str_list = [f"{int(ni)}/{int(N)}" for ni in n]
    
    if q == 1:
        # Special case for q=1: D_1 = exp(-sum(p_i * ln(p_i)))
        # This is the exponential of the Shannon entropy
        h_sum_terms = []
        for pi in p:
            if pi > 0:
                h_sum_terms.append(pi * np.log(pi))
        shannon_entropy = -np.sum(h_sum_terms)
        result = np.exp(shannon_entropy)
        
        # Format the equation string to be printed
        equation_str = "exp( - ( "
        log_terms_str = []
        for prop_str in proportions_str_list:
            log_terms_str.append(f"({prop_str}) * ln({prop_str})")
        equation_str += " + ".join(log_terms_str)
        equation_str += " ) )"
        print(f"Formula: D_1 = {equation_str}")

    else:
        # General formula for q != 1: D_q = (sum(p_i^q))^(1/(1-q))
        sum_p_q = np.sum(p**q)
        # Avoid division by zero if sum_p_q is zero
        if sum_p_q == 0:
            result = 0.0
        else:
            exponent = 1.0 / (1.0 - q)
            result = sum_p_q**exponent
        
        # Format the equation string to be printed
        equation_str = "( "
        power_terms_str = []
        for prop_str in proportions_str_list:
            power_terms_str.append(f"({prop_str})^{q}")
        equation_str += " + ".join(power_terms_str)
        equation_str += f" ) ^ (1 / (1 - {q}))"
        print(f"Formula: D_{q} = {equation_str}")
    
    print(f"Result for {level_name}: {result:.2f}")
    print("-" * (len(level_name) + 26))
    return result

# --- Step 1: Define abundances based on visual identification ---
# Based on visual inspection, we identified 4 distinct morphospecies with the following counts:
species_counts = [75, 15, 11, 6]

# --- Step 2: Define abundances at each taxonomic level based on the plan ---
# Order Level (all Lepidoptera): All individuals in one group.
order_abundances = [sum(species_counts)]  # [107]
# Family Level (all Nymphalidae): All individuals in one group.
family_abundances = [sum(species_counts)] # [107]
# Genus Level: 3 species in one genus, 1 in another.
# Genus 1 (likely Adelpha) = 75 + 15 + 11 = 101. Genus 2 = 6.
genus_abundances = [101, 6]
# Species Level: The 4 identified morphospecies.
species_abundances = species_counts

# --- Step 3: Calculate Hill numbers for each required case ---
hill_order = calculate_hill_number(order_abundances, 1, "Order")
hill_family = calculate_hill_number(family_abundances, 2, "Family")
hill_genus = calculate_hill_number(genus_abundances, 3, "Genus")
hill_species = calculate_hill_number(species_abundances, 4, "Species")

# --- Step 4: Print the final combined result ---
final_results = [hill_order, hill_family, hill_genus, hill_species]
formatted_results = ", ".join([f"{r:.2f}" for r in final_results])
print("\nFinal comma-separated values:")
print(formatted_results)

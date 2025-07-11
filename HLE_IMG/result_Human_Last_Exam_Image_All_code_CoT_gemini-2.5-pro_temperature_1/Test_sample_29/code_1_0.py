import numpy as np

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates the Hill number for a given list of abundances and order q.
    """
    # Convert to a numpy array for vectorized operations
    abundances = np.array(abundances, dtype=float)
    
    # Calculate total number of individuals
    total_individuals = np.sum(abundances)
    
    # Calculate proportions (p_i)
    proportions = abundances / total_individuals
    
    # Filter out zero proportions to avoid errors in log
    proportions = proportions[proportions > 0]
    
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    print(f"Abundances (n_i): {abundances.astype(int)}")
    print(f"Proportions (p_i): {[round(p, 4) for p in proportions]}")
    
    if q == 1:
        # Special case for q=1: exp(Shannon entropy)
        shannon_entropy_sum = -np.sum(proportions * np.log(proportions))
        hill_val = np.exp(shannon_entropy_sum)
        print(f"Equation: exp(-Σ(p_i * ln(p_i))) = exp({round(shannon_entropy_sum, 4)})")
    else:
        # General formula for q != 1
        sum_of_powers = np.sum(proportions**q)
        exponent = 1 / (1 - q)
        hill_val = sum_of_powers**exponent
        print(f"Equation: (Σ(p_i^{q}))^(1/(1-q)) = ({round(sum_of_powers, 4)})^({round(exponent, 4)})")
        
    print(f"Result for {level_name}: {round(hill_val, 2)}\n")
    return hill_val

# --- Data based on image analysis ---

# 1. Order Level: All 86 specimens belong to one order (Lepidoptera)
abundances_order = [86]
hill_order = calculate_hill_number(abundances_order, q=1, level_name="Order")

# 2. Family Level: All 86 specimens belong to one family (Nymphalidae)
abundances_family = [86]
hill_family = calculate_hill_number(abundances_family, q=2, level_name="Family")

# 3. Genus Level: Two genera identified with counts 78 and 8
abundances_genus = [78, 8]
hill_genus = calculate_hill_number(abundances_genus, q=3, level_name="Genus")

# 4. Species Level: Four species identified with counts 54, 6, 18, and 8
abundances_species = [54, 6, 18, 8]
hill_species = calculate_hill_number(abundances_species, q=4, level_name="Species")

# --- Final Output ---
final_results = [hill_order, hill_family, hill_genus, hill_species]
formatted_results = ", ".join([f"{res:.2f}" for res in final_results])

print("Final comma-separated values:")
print(formatted_results)

# Final answer in the specified format
print(f"<<<{formatted_results}>>>")
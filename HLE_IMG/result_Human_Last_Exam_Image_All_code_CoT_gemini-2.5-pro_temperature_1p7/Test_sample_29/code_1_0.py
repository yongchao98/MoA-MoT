import math

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates and prints the Hill number of order q for a given list of abundances.
    """
    total_individuals = sum(abundances)
    proportions = [n / total_individuals for n in abundances]

    # Print the "equation" components
    print(f"Calculating Hill Number for {level_name} level with q={q}:")
    print(f"Abundance counts (n_i): {abundances}")
    print(f"Proportions (p_i): {[round(p, 4) for p in proportions]}")

    if q == 1:
        # Special case for q=1: exponential of Shannon entropy
        if any(p == 1.0 for p in proportions):
             # The sum is 1*ln(1) = 0, so D_1 = exp(0) = 1
            result = 1.0
        else:
            shannon_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
            result = math.exp(shannon_entropy)
    else:
        # General formula for q != 1
        sum_p_q = sum(p**q for p in proportions if p > 0)
        result = sum_p_q**(1 / (1 - q))

    print(f"Result (D_{q}): {result:.4f}\n")
    return result

# --- Step 1: Define abundance counts based on image analysis ---

# Counts for the 3 identified species
species_1_count = 90
species_2_count = 19
species_3_count = 6
total_count = species_1_count + species_2_count + species_3_count

# Define abundance lists for each taxonomic level
abundances_order = [total_count]
abundances_family = [total_count]
abundances_genus = [species_1_count + species_2_count, species_3_count]
abundances_species = [species_1_count, species_2_count, species_3_count]

# --- Step 2: Calculate the four Hill numbers ---

hill_order_q1 = calculate_hill_number(abundances_order, 1, "Order")
hill_family_q2 = calculate_hill_number(abundances_family, 2, "Family")
hill_genus_q3 = calculate_hill_number(abundances_genus, 3, "Genus")
hill_species_q4 = calculate_hill_number(abundances_species, 4, "Species")

# --- Step 3: Format and print the final answer ---

final_answer_string = (
    f"{hill_order_q1:.2f}, {hill_family_q2:.2f}, "
    f"{hill_genus_q3:.2f}, {hill_species_q4:.2f}"
)

print(f"Final comma-separated values: {final_answer_string}")
print(f"<<<{final_answer_string}>>>")
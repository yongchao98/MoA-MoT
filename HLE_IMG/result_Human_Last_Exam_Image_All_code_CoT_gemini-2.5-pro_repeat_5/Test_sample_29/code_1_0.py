import numpy as np

def calculate_and_print_hill(level_name, counts, q):
    """
    Calculates the Hill number and prints the steps of the calculation.
    Returns the numeric result.
    """
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    
    counts = np.array(counts, dtype=float)
    total_individuals = np.sum(counts)
    proportions = counts / total_individuals
    
    print(f"Counts for each group: {counts.astype(int).tolist()}")
    print(f"Total individuals (N): {int(total_individuals)}")
    print(f"Proportions (p_i): {[round(p, 4) for p in proportions]}")
    
    # Handle the case of only one group, where diversity is always 1
    if len(counts) <= 1:
        result = 1.0
        print("Since there is only one group, the Hill number is 1.00\n")
        return result

    if q == 1:
        # Hill number for q=1 is the exponential of Shannon entropy
        shannon_entropy = -np.sum(proportions * np.log(proportions))
        result = np.exp(shannon_entropy)
        print(f"Shannon Entropy (-Σ p_i * ln(p_i)): {shannon_entropy:.4f}")
        print(f"Final Hill Number (exp(Entropy)): {result:.2f}\n")
    else:
        # General formula for q != 1
        sum_p_q = np.sum(proportions ** q)
        exponent = 1 / (1 - q)
        result = sum_p_q ** exponent
        print(f"Sum of proportions to the power of q (Σ p_i^{q}): {sum_p_q:.4f}")
        print(f"Final Hill Number ((Σ p_i^{q}) ^ (1/(1-q))): {result:.2f}\n")
        
    return result

# Based on visual inspection of the image, the counts for each taxonomic level are:
order_counts = [96]
family_counts = [96]
genus_counts = [84, 6, 6]
species_counts = [48, 16, 20, 6, 6]

# Calculate the four Hill numbers as requested
h_order = calculate_and_print_hill("Order", order_counts, 1)
h_family = calculate_and_print_hill("Family", family_counts, 2)
h_genus = calculate_and_print_hill("Genus", genus_counts, 3)
h_species = calculate_and_print_hill("Species", species_counts, 4)

# Print the final answer in the required format
final_answer = f"{h_order:.2f},{h_family:.2f},{h_genus:.2f},{h_species:.2f}"
print("--- Final Result ---")
print("The Hill numbers for order (q=1), family (q=2), genus (q=3), and species (q=4) are:")
print(final_answer)
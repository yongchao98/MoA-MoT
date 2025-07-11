import numpy as np

def calculate_hill_number(counts, q):
    """
    Calculates the Hill number of order q for a given set of counts.

    Args:
        counts (list of int): A list of abundances for each group.
        q (float or int): The order of the Hill number.

    Returns:
        float: The calculated Hill number.
    """
    # Convert counts to proportions
    total_individuals = sum(counts)
    if total_individuals == 0:
        return 0.0
    
    proportions = np.array([count / total_individuals for count in counts])
    
    # Filter out zero proportions to avoid math errors (e.g., log(0))
    proportions = proportions[proportions > 0]
    
    if q == 1:
        # Special case for q=1 (Exponential of Shannon entropy)
        shannon_entropy = -np.sum(proportions * np.log(proportions))
        return np.exp(shannon_entropy)
    else:
        # General formula for q != 1
        sum_p_q = np.sum(proportions ** q)
        if sum_p_q == 0:
             return 0.0
        return sum_p_q ** (1 / (1 - q))

def main():
    """
    Main function to calculate and print Hill numbers for the insect collection.
    """
    # Counts based on visual analysis of the image
    # Total individuals = 65 + 20 + 5 + 1 = 91
    counts_order = [91]           # All belong to one Order (Lepidoptera)
    counts_family = [91]          # All belong to one Family (Nymphalidae)
    counts_genus = [65, 20, 6]      # Three identified genera
    counts_species = [65, 20, 5, 1] # Four identified species

    # q values for each taxonomic level
    q_order = 1
    q_family = 2
    q_genus = 3
    q_species = 4
    
    # Calculate Hill numbers
    hill_order = calculate_hill_number(counts_order, q_order)
    hill_family = calculate_hill_number(counts_family, q_family)
    hill_genus = calculate_hill_number(counts_genus, q_genus)
    hill_species = calculate_hill_number(counts_species, q_species)
    
    # Print the results as requested
    results = [hill_order, hill_family, hill_genus, hill_species]
    print(f"The Hill number for Order (q=1) is: {hill_order:.2f}")
    print(f"The Hill number for Family (q=2) is: {hill_family:.2f}")
    print(f"The Hill number for Genus (q=3) is: {hill_genus:.2f}")
    print(f"The Hill number for Species (q=4) is: {hill_species:.2f}")
    print("\nFinal comma-separated values:")
    print(','.join([f'{res:.2f}' for res in results]))

if __name__ == "__main__":
    main()
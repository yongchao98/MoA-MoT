import numpy as np

def calculate_hill_numbers():
    """
    This script calculates Hill numbers for a sample of butterflies based on visual identification from an image.
    """

    # Step 1: Data from visual inspection of the image.
    # The specimens are classified into four distinct morphospecies (M1, M2, M3, M4).
    species_counts = {
        'M1': 67,  # Dominant orange and black striped butterfly
        'M2': 22,  # Black butterfly with bold white bands
        'M3': 3,   # Dark butterfly with a complex orange/white pattern
        'M4': 1    # Unique individual with a broken orange band
    }
    species_list_counts = list(species_counts.values())
    total_individuals = sum(species_list_counts)

    # Step 2: Define taxonomic hierarchy and counts for each level.
    # Order Level: All butterflies are in the Order Lepidoptera.
    order_counts = [total_individuals]
    # Family Level: All are assumed to be in the Family Nymphalidae.
    family_counts = [total_individuals]
    # Genus Level: Assume two genera based on distinct visual patterns.
    # Genus 1 (M1, M3, M4) and Genus 2 (M2).
    genus_counts = [species_counts['M1'] + species_counts['M3'] + species_counts['M4'], species_counts['M2']]
    # Species Level: The four identified morphospecies.
    # The counts are already in species_list_counts.

    print("Calculating Hill numbers based on the following counts:")
    print(f"Total Individuals (N): {total_individuals}")
    print(f"Species Counts: {species_list_counts}")
    print(f"Genus Counts: {genus_counts}")
    print(f"Family Counts: {family_counts}")
    print(f"Order Counts: {order_counts}\n")


    # Step 3: Define Hill number calculation function.
    def get_hill_number(counts, q):
        """Calculates the Hill number of order q from a list of counts."""
        if sum(counts) == 0:
            return 0.0
        proportions = np.array(counts) / sum(counts)
        # Filter out zero proportions for log calculation
        proportions = proportions[proportions > 0]
        
        if q == 1:
            # Special case for q=1 (exponential of Shannon entropy)
            return np.exp(-np.sum(proportions * np.log(proportions)))
        else:
            # General formula for q != 1
            return np.sum(proportions ** q) ** (1 / (1 - q))

    # Step 4: Calculate Hill numbers for each required level and q value.
    # D(q=1) for Order
    d1_order = get_hill_number(order_counts, 1)
    p_order = order_counts[0] / total_individuals
    print("Hill Number for Order (q=1):")
    print(f"D(q=1) = exp( - (({order_counts[0]}/{total_individuals})*ln({order_counts[0]}/{total_individuals})) ) = {d1_order:.2f}\n")

    # D(q=2) for Family
    d2_family = get_hill_number(family_counts, 2)
    p_family = family_counts[0] / total_individuals
    print("Hill Number for Family (q=2):")
    print(f"D(q=2) = ( ({family_counts[0]}/{total_individuals})^2 ) ^ (1/(1-2)) = {d2_family:.2f}\n")

    # D(q=3) for Genus
    d3_genus = get_hill_number(genus_counts, 3)
    g1, g2 = genus_counts[0], genus_counts[1]
    print("Hill Number for Genus (q=3):")
    print(f"D(q=3) = ( ({g1}/{total_individuals})^3 + ({g2}/{total_individuals})^3 ) ^ (1/(1-3)) = {d3_genus:.2f}\n")

    # D(q=4) for Species
    d4_species = get_hill_number(species_list_counts, 4)
    s1, s2, s3, s4 = species_list_counts[0], species_list_counts[1], species_list_counts[2], species_list_counts[3]
    print("Hill Number for Species (q=4):")
    print(f"D(q=4) = ( ({s1}/{total_individuals})^4 + ({s2}/{total_individuals})^4 + ({s3}/{total_individuals})^4 + ({s4}/{total_individuals})^4 ) ^ (1/(1-4)) = {d4_species:.2f}\n")
    
    # Step 5: Format and print the final answer.
    final_answer = f"{d1_order:.2f},{d2_family:.2f},{d3_genus:.2f},{d4_species:.2f}"
    print("Final comma-separated values rounded to two decimal places:")
    print(final_answer)
    return final_answer

final_answer_string = calculate_hill_numbers()
print(f"<<<{final_answer_string}>>>")
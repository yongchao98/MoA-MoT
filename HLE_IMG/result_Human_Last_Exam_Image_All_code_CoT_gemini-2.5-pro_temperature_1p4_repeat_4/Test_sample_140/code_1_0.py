import math

def calculate_shannon_diversity(species_counts):
    """
    Calculates the Shannon's diversity index (H).
    H = -SUM(pi * ln(pi))
    """
    total_individuals = sum(species_counts)
    if total_individuals == 0:
        return 0.0

    shannon_index = 0.0
    for count in species_counts:
        # A site with only one species will have a proportion of 1, and ln(1) is 0, so H=0.
        if count > 0:
            proportion = count / total_individuals
            shannon_index += proportion * math.log(proportion)

    return -shannon_index

# Step 1 & 2: Define initial data and consolidate based on true species identity
# Morphospecies A, B, C, D, E are Harmonia axyridis
# Morphospecies F is a Tortoise Beetle

# Site 1
# Original: A: 1, B: 1, C: 2
site1_true_counts = [1 + 1 + 2] # All are H. axyridis

# Site 2
# Original: A: 10, D: 12, E: 1
site2_true_counts = [10 + 12 + 1] # All are H. axyridis

# Site 3
# Original: B: 2, E: 2, F: 4
site3_true_counts = [2 + 2, 4] # H. axyridis and Tortoise Beetle

# Step 3: Calculate diversity for each site
h1 = calculate_shannon_diversity(site1_true_counts)
h2 = calculate_shannon_diversity(site2_true_counts)
h3 = calculate_shannon_diversity(site3_true_counts)

# Step 4: Display the results
print("Based on the true species identities, the abundances at each site are:")
print(f"Site 1: Harmonia axyridis: {site1_true_counts[0]}")
print(f"Site 2: Harmonia axyridis: {site2_true_counts[0]}")
print(f"Site 3: Harmonia axyridis: {site3_true_counts[0]}, Tortoise Beetle: {site3_true_counts[1]}")
print("\nShannon's Diversity Calculation (H = -SUM(pi * ln(pi))):")
print("Site 1: Contains only one species, so H = 0.00")
print("Site 2: Contains only one species, so H = 0.00")
n1_site3 = site3_true_counts[0]
n2_site3 = site3_true_counts[1]
N_site3 = sum(site3_true_counts)
p1_site3 = n1_site3 / N_site3
p2_site3 = n2_site3 / N_site3
print(f"Site 3: H = -[({n1_site3}/{N_site3} * ln({n1_site3}/{N_site3})) + ({n2_site3}/{N_site3} * ln({n2_site3}/{N_site3}))] = {h3:.2f}")


final_answer = f"{h1:.2f}, {h2:.2f}, {h3:.2f}"
print(f"\nFinal Answer: {final_answer}")
print(f"<<<{final_answer}>>>")
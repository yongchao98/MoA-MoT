import numpy as np

def hill_number(abundances, q):
    """
    Calculates the Hill number (effective number of species) for a given 
    list of abundances and an order q.
    """
    # Convert list to a numpy array for efficient calculation
    abundances = np.array(abundances, dtype=float)
    
    # Calculate the total number of individuals
    total_individuals = np.sum(abundances)
    
    # If there are no individuals, diversity is 0
    if total_individuals == 0:
        return 0.0
        
    # Calculate the proportional abundance of each group
    proportions = abundances / total_individuals
    
    # Filter out groups with zero abundance to avoid errors in log(0)
    proportions = proportions[proportions > 0]
    
    # The formula for q=1 is the limit as q->1, which is exp(Shannon entropy)
    if q == 1:
        shannon_entropy = -np.sum(proportions * np.log(proportions))
        return np.exp(shannon_entropy)
    # The general formula for Hill numbers for all other q
    else:
        sum_of_proportions_to_q = np.sum(proportions**q)
        return sum_of_proportions_to_q**(1 / (1 - q))

# --- Step 1: Define abundances based on taxonomic hierarchy ---

# Order level (q=1): All 108 individuals belong to one order (Lepidoptera).
abundances_order = [108]
q_order = 1

# Family level (q=2): All 108 individuals belong to one family (Nymphalidae).
abundances_family = [108]
q_family = 2

# Genus level (q=3): Assuming 4 distinct species represent 4 distinct genera.
# The abundances of the four genera are:
abundances_genus = [78, 4, 9, 17]
q_genus = 3

# Species level (q=4): We identified 4 distinct species.
# The abundances of the four species are:
abundances_species = [78, 4, 9, 17]
q_species = 4

# --- Step 2: Calculate the Hill number for each case ---

# D(q=1) for the Order level
d1_order = hill_number(abundances_order, q_order)

# D(q=2) for the Family level
d2_family = hill_number(abundances_family, q_family)

# D(q=3) for the Genus level
d3_genus = hill_number(abundances_genus, q_genus)

# D(q=4) for the Species level
d4_species = hill_number(abundances_species, q_species)

# --- Step 3: Print the final answer ---
# The final result consists of the four calculated Hill numbers,
# rounded to two decimal places and separated by commas.
print(f"{d1_order:.2f},{d2_family:.2f},{d3_genus:.2f},{d4_species:.2f}")
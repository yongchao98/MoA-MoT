import math

def calculate_hill_number(counts, q):
    """Calculates the Hill number for a given list of abundance counts and order q."""
    
    # Calculate total number of individuals
    total_individuals = sum(counts)
    if total_individuals == 0:
        return 0.0

    # Calculate proportions
    proportions = [count / total_individuals for count in counts]

    # Handle the special case for q=1 (Shannon diversity)
    if q == 1:
        shannon_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
        return math.exp(shannon_entropy)
    
    # General formula for q != 1
    else:
        sum_of_proportions_raised_to_q = sum(p**q for p in proportions)
        if sum_of_proportions_raised_to_q == 0:
            return 0.0
        return sum_of_proportions_raised_to_q**(1 / (1 - q))

def solve_and_print():
    """
    Solves the Hill number problem based on the provided image and prints the results.
    """
    print("--- Analysis of Insect Collection ---")
    
    # Step 1 & 2: Counts and Taxonomic Hierarchy based on visual analysis.
    # Total count = 56+20+4+18+4+3+1+4 = 110 individuals.
    species_counts = [56, 20, 4, 18, 4, 3, 1, 4]
    
    # Groups are combined for higher taxonomic levels.
    # Family/Genus Group 1 (red/orange): 56+20+4 = 80
    # Family/Genus Group 2 (dark/banded): 18+4+3+1+4 = 30
    family_genus_counts = [80, 30]
    
    # All individuals belong to one order.
    order_counts = [sum(species_counts)]
    
    results = []

    # --- Calculation for Order (q=1) ---
    q_order = 1
    d_order = calculate_hill_number(order_counts, q_order)
    results.append(d_order)
    print(f"\n1. Taxon: Order, q = {q_order}")
    print(f"   Counts: {order_counts}")
    print(f"   Equation for q=1: exp(-sum(p_i * ln(p_i)))")
    print(f"   Resulting Hill Number D({q_order}): {d_order:.2f}")

    # --- Calculation for Family (q=2) ---
    q_family = 2
    d_family = calculate_hill_number(family_genus_counts, q_family)
    results.append(d_family)
    print(f"\n2. Taxon: Family, q = {q_family}")
    print(f"   Counts: {family_genus_counts}")
    print(f"   Equation: (sum(p_i^{q_family}))^(1/(1-{q_family}))")
    print(f"   Resulting Hill Number D({q_family}): {d_family:.2f}")

    # --- Calculation for Genus (q=3) ---
    q_genus = 3
    d_genus = calculate_hill_number(family_genus_counts, q_genus)
    results.append(d_genus)
    print(f"\n3. Taxon: Genus, q = {q_genus}")
    print(f"   Counts: {family_genus_counts}")
    print(f"   Equation: (sum(p_i^{q_genus}))^(1/(1-{q_genus}))")
    print(f"   Resulting Hill Number D({q_genus}): {d_genus:.2f}")

    # --- Calculation for Species (q=4) ---
    q_species = 4
    d_species = calculate_hill_number(species_counts, q_species)
    results.append(d_species)
    print(f"\n4. Taxon: Species, q = {q_species}")
    print(f"   Counts: {species_counts}")
    print(f"   Equation: (sum(p_i^{q_species}))^(1/(1-{q_species}))")
    print(f"   Resulting Hill Number D({q_species}): {d_species:.2f}")

    # --- Final formatted output ---
    final_answer = ",".join([f"{r:.2f}" for r in results])
    print("\n--- Final Answer ---")
    print("The four numeric values (for Order, Family, Genus, Species) rounded to two decimal places are:")
    print(final_answer)

solve_and_print()
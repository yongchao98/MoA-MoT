import math
import numpy as np

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates the Hill number for a given set of abundances and order q.
    Also prints the equation used for the calculation.
    """
    N = sum(abundances)
    if N == 0:
        return 0.0

    proportions = [count / N for count in abundances]
    
    print(f"Calculating Hill Number for {level_name} (q={q}):")
    
    if q == 1:
        # Special case for q=1: exp(Shannon entropy)
        entropy_terms_str = []
        shannon_entropy = 0
        for p in proportions:
            if p > 0:
                shannon_entropy -= p * math.log(p)
                entropy_terms_str.append(f"{p:.3f} * ln({p:.3f})")

        print(f"  Formula: D^1 = exp(-Σ(pi * ln(pi)))")
        print(f"  Equation: D^1 = exp(-({" + ".join(entropy_terms_str)}))")
        result = math.exp(shannon_entropy)
        print(f"  Result: {result:.4f}\n")
        return result
    else:
        # General case for q != 1
        sum_p_q = sum(p**q for p in proportions)
        
        # Build the equation string
        equation_parts = [f"({count}/{N})^{q}" for count in abundances]
        
        print(f"  Formula: D^{q} = (Σ(pi^{q}))^(1/(1-{q}))")
        print(f"  Equation: D^{q} = ({" + ".join(equation_parts)})^(1/{1-q})")

        # Handle potential floating point issues if sum_p_q is zero or negative
        if sum_p_q <= 0:
            result = 0.0
        else:
            result = sum_p_q**(1 / (1 - q))

        print(f"  Result: {result:.4f}\n")
        return result

# Abundance data based on visual inspection of the butterfly collection
order_abundances = [113]
family_abundances = [113]
genus_abundances = [107, 6]
species_abundances = [82, 21, 4, 4, 1, 1]

# q values as specified in the prompt
q_order = 1
q_family = 2
q_genus = 3
q_species = 4

# Calculate Hill numbers for each level
h1_order = calculate_hill_number(order_abundances, q_order, "Order")
h2_family = calculate_hill_number(family_abundances, q_family, "Family")
h3_genus = calculate_hill_number(genus_abundances, q_genus, "Genus")
h4_species = calculate_hill_number(species_abundances, q_species, "Species")

# Format and print the final answer
final_answer = f"{h1_order:.2f},{h2_family:.2f},{h3_genus:.2f},{h4_species:.2f}"
print("Final comma-separated values (rounded to two decimal places):")
print(final_answer)
print(f"<<<{final_answer}>>>")
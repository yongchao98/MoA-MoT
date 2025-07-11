import math

def calculate_hill_number(counts, q, level_name):
    """Calculates the Hill number for a given set of counts and q value."""
    
    total = sum(counts)
    if total == 0:
        return 0.0

    proportions = [count / total for count in counts]
    
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    print(f"Counts: {counts}, Total Individuals: {total}")
    
    if q == 1:
        # Special case for q=1: exponential of Shannon entropy
        shannon_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
        equation_str_parts = []
        for c in counts:
            if c > 0:
                equation_str_parts.append(f"({c}/{total})*ln({c}/{total})")
        equation_str = f"exp( - ( " + " + ".join(equation_str_parts) + " ) )"
        print(f"Formula: ^1D = {equation_str}")
        result = math.exp(shannon_entropy)
    else:
        # General formula for q != 1
        sum_of_powers = sum(p**q for p in proportions)
        equation_str_parts = []
        for c in counts:
            equation_str_parts.append(f"({c}/{total})^{q}")
        equation_str = f"( " + " + ".join(equation_str_parts) + f" )^(1/(1-{q}))"
        print(f"Formula: ^{q}D = {equation_str}")
        
        # Avoid division by zero if sum_of_powers is 0
        if sum_of_powers == 0:
            result = 0.0
        else:
            result = sum_of_powers**(1 / (1 - q))
            
    print(f"Result for {level_name}: {result:.2f}\n")
    return result

# Data from image analysis
order_counts = [93]
family_counts = [93]
# Genus counts: Species A+B (66), Species C (21), Species D (6)
genus_counts = [66, 21, 6] 
# Species counts: Species A (54), C (21), B (12), D (6)
species_counts = [54, 21, 12, 6] 

# q values as per the problem
q_values = {'Order': 1, 'Family': 2, 'Genus': 3, 'Species': 4}

# Calculate Hill numbers
h1_order = calculate_hill_number(order_counts, q_values['Order'], 'Order')
h2_family = calculate_hill_number(family_counts, q_values['Family'], 'Family')
h3_genus = calculate_hill_number(genus_counts, q_values['Genus'], 'Genus')
h4_species = calculate_hill_number(species_counts, q_values['Species'], 'Species')

# Format final output as requested
final_answer = f"{h1_order:.2f},{h2_family:.2f},{h3_genus:.2f},{h4_species:.2f}"
print("--- Final Answer ---")
print("The Hill numbers for order (q=1), family (q=2), genus (q=3), and species (q=4) are:")
print(final_answer)

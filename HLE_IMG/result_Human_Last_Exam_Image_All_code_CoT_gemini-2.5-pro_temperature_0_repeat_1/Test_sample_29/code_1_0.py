import math

def calculate_hill_number(counts, q, level_name):
    """
    Calculates the Hill number for a given set of counts and order q,
    and prints the detailed calculation steps.
    """
    N = sum(counts)
    if N == 0:
        print(f"No individuals for {level_name}, Hill number is 0.")
        return 0.0
    
    proportions = [n / N for n in counts]
    
    print(f"Calculating Hill Number for {level_name} (q={q})")
    print(f"  - Counts: {counts}")
    print(f"  - Total (N): {N}")
    print(f"  - Proportions (p_i): {[round(p, 4) for p in proportions]}")

    if q == 1:
        # Special case for q=1 (Shannon diversity)
        print(f"  - Formula: D^1 = exp(-Σ(p_i * ln(p_i)))")
        # For a single group, p_1=1, and 1*ln(1)=0.
        shannon_entropy = -sum(p * math.log(p) for p in proportions if p > 0)
        result = math.exp(shannon_entropy)
        print(f"  - Calculation: exp(-({proportions[0]:.4f} * ln({proportions[0]:.4f}))) = exp(0)")
        print(f"  - Result: {result:.2f}\n")
        return result
    else:
        # General case for q != 1
        print(f"  - Formula: D^q = (Σ(p_i^q))^(1 / (1 - q))")
        sum_p_q = sum(p**q for p in proportions)
        
        # Build the calculation string
        calc_str_parts = [f"{p:.4f}^{q}" for p in proportions]
        calc_str = " + ".join(calc_str_parts)
        exponent = 1 / (1 - q)
        
        print(f"  - Calculation: ({calc_str}) ^ (1 / (1 - {q}))")
        print(f"               = ({sum_p_q:.4f}) ^ ({exponent:.4f})")
        
        if sum_p_q == 0:
            result = 0.0
        else:
            result = sum_p_q**exponent
        
        print(f"  - Result: {result:.2f}\n")
        return result

# Step 1: Define counts based on visual classification.
counts_species = [62, 4, 18, 7] # Species A, B, C, D
total_individuals = sum(counts_species)

# Step 2: Define counts for each taxonomic level.
counts_order = [total_individuals]
counts_family = [62 + 4 + 7, 18] # [73, 18]
counts_genus = [62 + 4, 18, 7]   # [66, 18, 7]

# Step 3: Define q values.
q_order = 1
q_family = 2
q_genus = 3
q_species = 4

# Step 4: Calculate Hill numbers for each level.
hill_order = calculate_hill_number(counts_order, q_order, "Order")
hill_family = calculate_hill_number(counts_family, q_family, "Family")
hill_genus = calculate_hill_number(counts_genus, q_genus, "Genus")
hill_species = calculate_hill_number(counts_species, q_species, "Species")

# Step 5: Format and print the final answer.
final_answer = f"{hill_order:.2f},{hill_family:.2f},{hill_genus:.2f},{hill_species:.2f}"
print("Final comma-separated values:")
print(final_answer)
import math

def calculate_shannon_diversity(site_name, counts):
    """Calculates Shannon's diversity and prints the equation and result."""
    total_individuals = sum(counts)
    
    if total_individuals == 0 or len(counts) <= 1:
        shannon_index = 0.0
        # Build equation string for a single species
        if total_individuals > 0:
            count = counts[0]
            equation_str = f"H' = - [ ({count}/{total_individuals}) * ln({count}/{total_individuals}) ]"
        else: # Should not happen with the given data
            equation_str = "H' = 0 (No individuals)"
    else:
        shannon_sum = 0
        equation_parts = []
        for count in counts:
            if count > 0:
                proportion = count / total_individuals
                shannon_sum += proportion * math.log(proportion)
                equation_parts.append(f"({count}/{total_individuals}) * ln({count}/{total_individuals})")
        
        shannon_index = -shannon_sum
        equation_str = f"H' = - [ {' + '.join(equation_parts)} ]"

    print(f"{site_name}:")
    print(f"Abundances: {counts}")
    print(f"Calculation: {equation_str} = {shannon_index:.2f}")
    return shannon_index

# Corrected insect abundance data for each site
site1_abundances = [4]
site2_abundances = [11, 12]
site3_abundances = [4]

# Calculate and print diversity for each site
h1 = calculate_shannon_diversity("Site 1", site1_abundances)
print("-" * 20)
h2 = calculate_shannon_diversity("Site 2", site2_abundances)
print("-" * 20)
h3 = calculate_shannon_diversity("Site 3", site3_abundances)
print("-" * 20)

# Final answer formatting
final_answer = f"{h1:.2f}, {h2:.2f}, {h3:.2f}"
print(f"\nThe insect diversity of each site is: {final_answer}")

# Final result in the specified format
print(f"<<<{final_answer}>>>")
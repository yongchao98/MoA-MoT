import math

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates the Hill number for a given set of abundances and order q.
    It also prints the detailed calculation steps.
    """
    total_individuals = sum(abundances)
    if total_individuals == 0:
        return 0.0

    proportions = [n / total_individuals for n in abundances]
    
    print(f"Calculating Hill number for {level_name} with q={q}")
    print(f"Abundances: {abundances}, Total Individuals: {total_individuals}")

    if q == 1:
        # Special case for q=1: exponential of Shannon entropy
        shannon_entropy_sum_str = []
        shannon_entropy_sum = 0
        for p in proportions:
            if p > 0:
                shannon_entropy_sum -= p * math.log(p)
                shannon_entropy_sum_str.append(f"{p:.2f}*ln({p:.2f})")
        
        equation_str = f"exp(-({' + '.join(shannon_entropy_sum_str)}))"
        result = math.exp(shannon_entropy_sum)
        print(f"Equation: {equation_str}")
        print(f"Result: {result:.2f}\n")
        return result

    else:
        # General case for q != 1
        sum_of_proportions_to_q = sum(p**q for p in proportions)
        
        # Build the equation string for printing
        sum_pq_str_parts = [f"({n}/{total_individuals})^{q}" for n in abundances]
        sum_pq_str = " + ".join(sum_pq_str_parts)
        exponent_str = f"1/(1-{q})"
        equation_str = f"({sum_pq_str}) ^ ({exponent_str})"
        
        if sum_of_proportions_to_q > 0:
            result = sum_of_proportions_to_q**(1 / (1 - q))
        else:
            result = 0.0
            
        print(f"Equation: {equation_str}")
        print(f"Result: {result:.2f}\n")
        return result

def main():
    """
    Main function to define data and calculate Hill numbers.
    """
    # Step 1 & 2: Data Extraction and Taxonomic Grouping
    # Counts are based on visual inspection of the insect trays.
    
    # Species abundances (5 morpho-species identified)
    # s1 (dominant orange-banded): 64
    # s2 (white-banded): 27 (11+11+5)
    # s3 (intricate pattern): 6
    # s4 (small orange-banded): 4
    # s5 (miscellaneous): 5
    species_abundances = [64, 27, 6, 4, 5] # Total = 106

    # Genus abundances (4 genera identified by grouping similar species)
    # Genus A (s1+s4): 64+4 = 68
    # Genus B (s2): 27
    # Genus C (s3): 6
    # Genus D (s5): 5
    genus_abundances = [68, 27, 6, 5] # Total = 106

    # Family abundances (All specimens appear to be in the same family)
    family_abundances = [106] # Total = 106
    
    # Order abundances (All specimens are Lepidoptera)
    order_abundances = [106] # Total = 106
    
    # Step 3 & 4: Calculate Hill numbers for each level with its respective q
    results = []
    
    # Order, q=1
    d_order = calculate_hill_number(order_abundances, 1, "Order")
    results.append(d_order)
    
    # Family, q=2
    d_family = calculate_hill_number(family_abundances, 2, "Family")
    results.append(d_family)

    # Genus, q=3
    d_genus = calculate_hill_number(genus_abundances, 3, "Genus")
    results.append(d_genus)

    # Species, q=4
    d_species = calculate_hill_number(species_abundances, 4, "Species")
    results.append(d_species)

    # Step 5: Format and print the final output
    final_output = ", ".join([f"{r:.2f}" for r in results])
    print("---")
    print("Final comma-separated values:")
    print(final_output)


if __name__ == "__main__":
    main()
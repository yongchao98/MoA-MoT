import numpy as np

def calculate_hill_number(abundances, q, level_name):
    """
    Calculates the Hill number for a given list of abundances and order q.
    """
    total_individuals = sum(abundances)
    proportions = np.array([n / total_individuals for n in abundances if n > 0])
    
    # Format the numbers used in the equation
    abundance_str = ", ".join(map(str, abundances))
    equation_str = f"For {level_name} level (q={q}) with abundances [{abundance_str}]:\n"
    
    if q == 1:
        # Handle the special case for q=1 (Shannon diversity)
        if len(proportions) == 1 and proportions[0] == 1.0:
            shannon_entropy = 0.0
        else:
            shannon_entropy = -np.sum(proportions * np.log(proportions))
        
        result = np.exp(shannon_entropy)
        
        pi_log_pi_terms = [f"{p:.2f}*ln({p:.2f})" for p in proportions]
        equation_str += f"  Formula: exp(- sum(pi * ln(pi)))\n"
        equation_str += f"  Calculation: exp(- ({' + '.join(pi_log_pi_terms)})) = {result:.2f}"

    else:
        # General formula for q != 1
        sum_of_powers = np.sum(proportions ** q)
        # To avoid division by zero if all proportions are zero (though filtered earlier)
        if sum_of_powers == 0:
            result = 0.0
        else:
            result = sum_of_powers ** (1 / (1 - q))
        
        pi_q_terms = [f"{p:.2f}^{q}" for p in proportions]
        equation_str += f"  Formula: (sum(pi^q))^(1/(1-q))\n"
        equation_str += f"  Calculation: ({' + '.join(pi_q_terms)})^(1/(1-{q})) = {result:.2f}"

    print(equation_str)
    return result

def solve():
    """
    Main function to solve the problem based on image analysis.
    """
    # Step 1 & 2: Counts and taxonomic grouping based on visual analysis
    # Species counts: A=87, B=21, C=6, D=1. Total=115
    species_abundances = [87, 21, 6, 1]

    # Genus grouping: Genus1(A), Genus2(B), Genus3(C+D)
    genus_abundances = [87, 21, 7]
    
    # Family grouping: Family1(A+C), Family2(B+D)
    family_abundances = [93, 22]

    # Order grouping: All individuals are in one order
    order_abundances = [115]

    # Step 3 & 4: Calculate Hill numbers for each level and q
    # q=1 for order, q=2 for family, q=3 for genus, q=4 for species
    d_order = calculate_hill_number(order_abundances, 1, "Order")
    d_family = calculate_hill_number(family_abundances, 2, "Family")
    d_genus = calculate_hill_number(genus_abundances, 3, "Genus")
    d_species = calculate_hill_number(species_abundances, 4, "Species")

    # Format the final output string
    final_answer = f"{d_order:.2f}, {d_family:.2f}, {d_genus:.2f}, {d_species:.2f}"
    
    print("\n--- Final Answer ---")
    print(final_answer)

solve()
<<<1.00, 1.45, 1.51, 1.45>>>
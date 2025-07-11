import numpy as np

def calculate_and_print_hill(level_name, abundances, q):
    """
    Calculates the Hill number for a given list of abundances and q value,
    and prints the steps of the calculation.
    """
    print(f"--- Calculating Hill Number for {level_name} (q={q}) ---")
    
    total_individuals = sum(abundances)
    
    # Handle case with only one group
    if len(abundances) == 1:
        print(f"Abundances: {abundances}")
        print("There is only one group, so the diversity is 1.")
        result = 1.0
        print(f"Result: {result:.2f}\n")
        return result

    proportions = np.array([n / total_individuals for n in abundances])
    proportions_str = ", ".join([f"{p:.2f}" for p in proportions])
    
    print(f"Abundances: {abundances}")
    print(f"Total Individuals (N): {total_individuals}")
    print(f"Proportions (p_i): [{proportions_str}]")

    # Handle q=1 case (Shannon diversity)
    if q == 1:
        proportions_gt_zero = proportions[proportions > 0]
        shannon_entropy = -np.sum(proportions_gt_zero * np.log(proportions_gt_zero))
        result = np.exp(shannon_entropy)
        print(f"Shannon Entropy (H): -sum(p_i * ln(p_i)) = {shannon_entropy:.4f}")
        print(f"Hill Number (D_1) = exp(H)")
    # General case for q != 1
    else:
        sum_p_q = np.sum(proportions ** q)
        result = sum_p_q ** (1 / (1 - q))
        print(f"Equation: (sum(p_i^{q}))^(1/(1-{q}))")
        print(f"Sum of p_i^{q}: {sum_p_q:.4f}")
        print(f"Hill Number (D_{q})")

    print(f"Result: {result:.2f}\n")
    return result

def solve_diversity():
    """
    Main function to solve the problem based on visual counts from the image.
    """
    # Step 1: Count individuals for each visually distinct species from the image.
    # Species A (orange/black, prominent bands): 12+12+8+6+12+12+6 = 68
    # Species B (dark/white stripes): 12+10 = 22
    # Species C (dark/orange, different pattern): 6
    count_A = 68
    count_B = 22
    count_C = 6

    # Step 2: Define abundance lists for each taxonomic level.
    # Order: All individuals in one group.
    abundances_order = [count_A + count_B + count_C]
    # Family: Two visually distinct groups (A+C vs B).
    abundances_family = [count_A + count_C, count_B]
    # Genus: Three visually distinct groups.
    abundances_genus = [count_A, count_B, count_C]
    # Species: Three visually distinct groups.
    abundances_species = [count_A, count_B, count_C]

    # Step 3: Calculate Hill numbers for each level with its specified q.
    hill_order = calculate_and_print_hill("Order", abundances_order, q=1)
    hill_family = calculate_and_print_hill("Family", abundances_family, q=2)
    hill_genus = calculate_and_print_hill("Genus", abundances_genus, q=3)
    hill_species = calculate_and_print_hill("Species", abundances_species, q=4)

    # Step 4: Format and print the final answer.
    final_results = [hill_order, hill_family, hill_genus, hill_species]
    final_string = ",".join([f"{res:.2f}" for res in final_results])
    
    print("--- Final Answer ---")
    print("The Hill numbers for order (q=1), family (q=2), genus (q=3), and species (q=4) are:")
    print(final_string)
    print(f"<<<{final_string}>>>")

# Run the solver
solve_diversity()
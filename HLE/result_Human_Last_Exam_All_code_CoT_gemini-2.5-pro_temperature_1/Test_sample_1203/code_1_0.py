import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which
    is most likely to lead to the rejection of the independent assortment hypothesis.
    """
    # Phenotype order:
    # 1. Tall, round, yellow (T_R_Y_)
    # 2. Tall, round, green (T_R_yy)
    # 3. Tall, wrinkled, yellow (T_rrY_)
    # 4. Tall, wrinkled, green (T_rryy)
    # 5. Dwarf, round, yellow (ttR_Y_)
    # 6. Dwarf, round, green (ttR_yy)
    # 7. Dwarf, wrinkled, yellow (ttrrY_)
    # 8. Dwarf, wrinkled, green (ttrryy)

    # Expected ratio for a trihybrid cross with independent assortment
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    ratio_total = np.sum(expected_ratio)

    # Observed counts for each answer choice
    observed_counts = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }
    
    phenotype_labels = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow", 
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green", 
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    results = {}

    print("Calculating Chi-Square values for each option:\n")

    for option, observed in observed_counts.items():
        total_offspring = np.sum(observed)
        expected_counts = (expected_ratio / ratio_total) * total_offspring
        
        # Chi-square test cannot have an expected value of 0.
        # If an expected value is very low, it's common to pool categories,
        # but for this problem, we'll proceed with the direct calculation.
        # Add a small epsilon to avoid division by zero if an expected count is 0.
        epsilon = 1e-9
        chi_square_terms = (observed - expected_counts)**2 / (expected_counts + epsilon)
        chi_square_value = np.sum(chi_square_terms)
        results[option] = chi_square_value
        print(f"Option {option}: Chi-Square = {chi_square_value:.2f}")

    # Find the option with the maximum chi-square value
    most_likely_option = max(results, key=results.get)
    max_chi_square = results[most_likely_option]

    print(f"\nOption {most_likely_option} has the highest Chi-Square value ({max_chi_square:.2f}), indicating the greatest deviation from the expected 27:9:9:3:9:3:3:1 ratio. This option would most likely lead to the rejection of the null hypothesis.")

    print(f"\nDetailed calculation for Option {most_likely_option}:")
    
    obs_final = observed_counts[most_likely_option]
    total_final = np.sum(obs_final)
    exp_final = (expected_ratio / ratio_total) * total_final
    
    equation_parts = []
    for i in range(len(obs_final)):
        o = obs_final[i]
        e = exp_final[i]
        term = (o - e)**2 / (e + epsilon)
        part_str = f"({o} - {e:.2f})^2 / {e:.2f}"
        equation_parts.append(part_str)
        print(f"  - For '{phenotype_labels[i]}': {part_str} = {term:.2f}")

    print("\nFinal Chi-Square Equation:")
    print(" + ".join(equation_parts) + f" = {max_chi_square:.2f}")
    
    print("\n<<<E>>>")


solve_chi_square()
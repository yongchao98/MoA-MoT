import sys

def solve_chi_square():
    """
    Calculates the chi-square statistic for each of the five answer choices
    to determine which is most likely to be rejected.
    """
    # Phenotypic classes in order:
    # 1. Tall, round, yellow (T_R_Y_)
    # 2. Tall, round, green (T_R_yy)
    # 3. Tall, wrinkled, yellow (T_rrY_)
    # 4. Tall, wrinkled, green (T_rryy)
    # 5. dwarf, round, yellow (ttR_Y_)
    # 6. dwarf, round, green (ttR_yy)
    # 7. dwarf, wrinkled, yellow (ttrrY_)
    # 8. dwarf, wrinkled, green (ttrryy)

    # Expected ratio for a trihybrid cross
    ratio_parts = [27, 9, 9, 9, 3, 3, 3, 1]
    total_ratio_parts = sum(ratio_parts)

    # Observed counts for each answer choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}
    
    # Calculate chi-square for each option
    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        
        # When total_observed is zero, chi-square cannot be calculated.
        if total_observed == 0:
            results[option] = (0, []) # Chi-square is 0, no expected values
            continue

        expected_counts = [(count / total_ratio_parts) * total_observed for count in ratio_parts]
        
        chi_square = 0
        for obs, exp in zip(observed_counts, expected_counts):
            # If expected is 0, the contribution is undefined. This shouldn't happen with this ratio.
            if exp == 0: 
                # To prevent division by zero, if obs is also 0, contribution is 0. If obs is not 0, it's infinite deviation.
                if obs > 0:
                    chi_square = float('inf')
                    break
                else:
                    continue

            chi_square += ((obs - exp)**2) / exp
            
        results[option] = (chi_square, expected_counts, observed_counts)

    # Print the results
    print("Chi-Square Test for Independent Assortment (Expected Ratio 27:9:9:9:3:3:3:1)\n")
    max_chi_square = -1
    best_option = None
    
    for option in sorted(results.keys()):
        chi_square_val, _, _ = results[option]
        print(f"Option {option}: χ² = {chi_square_val:.2f}")
        if chi_square_val > max_chi_square:
            max_chi_square = chi_square_val
            best_option = option

    print(f"\nConclusion: Option {best_option} has the highest chi-square value ({max_chi_square:.2f}).")
    print("This indicates the largest deviation from the expected ratio, making it the most likely to lead to a rejection of the null hypothesis of independent assortment.\n")

    # Display the full calculation for the winning option
    print(f"Detailed Calculation for Option {best_option}:")
    chi_square_val, expected_counts, observed_counts = results[best_option]
    
    equation_parts = []
    for i in range(len(observed_counts)):
        obs = observed_counts[i]
        exp = expected_counts[i]
        part_str = f"({obs} - {exp:.2f})² / {exp:.2f}"
        equation_parts.append(part_str)
        
    full_equation = " + ".join(equation_parts)
    print(f"χ² = {full_equation}")
    print(f"   = {chi_square_val:.2f}")

solve_chi_square()
print("<<<E>>>")
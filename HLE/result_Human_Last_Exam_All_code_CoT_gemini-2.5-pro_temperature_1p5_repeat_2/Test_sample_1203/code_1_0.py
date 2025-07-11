import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square statistic for five different sets of observed offspring
    phenotypes from a trihybrid cross to determine which set is most likely to
    lead to the rejection of the null hypothesis of independent assortment.
    """
    # Phenotype order:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. Dwarf, round, yellow
    # 6. Dwarf, round, green
    # 7. Dwarf, wrinkled, yellow
    # 8. Dwarf, wrinkled, green
    
    # Expected phenotypic ratio for a trihybrid cross
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)
    expected_ratio = [part / total_ratio_parts for part in expected_ratio_parts]

    # Observed data from answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}

    print("--- Chi-Square Calculation for Each Option ---")
    for choice, observed in observed_data.items():
        total_observed = sum(observed)
        if total_observed == 0:
            results[choice] = 0
            continue
            
        expected_counts = [ratio * total_observed for ratio in expected_ratio]
        
        # Calculate chi-square value
        chi_square_stat = sum(
            (o - e)**2 / e for o, e in zip(observed, expected_counts)
        )
        
        results[choice] = chi_square_stat
        print(f"Option {choice}: Total Offspring = {total_observed}, Chi-Square = {chi_square_stat:.2f}")
    
    # Find the choice with the maximum chi-square value
    best_choice = max(results, key=results.get)
    print(f"\nOption '{best_choice}' has the highest Chi-Square value, indicating the largest deviation from the expected ratio.")
    print("This option would most likely lead to the rejection of the null hypothesis.\n")
    
    # Print the detailed calculation for the best choice
    print(f"--- Detailed Calculation for Option {best_choice} ---")
    observed_final = observed_data[best_choice]
    total_observed_final = sum(observed_final)
    expected_final = [ratio * total_observed_final for ratio in expected_ratio]
    
    equation_parts = []
    for o, e in zip(observed_final, expected_final):
        equation_parts.append(f"({o} - {e:.2f})^2 / {e:.2f}")
        
    final_equation = " + ".join(equation_parts)
    print(f"Chi-Square = {final_equation}")
    print(f"Result = {results[best_choice]:.2f}\n")

    print("The degrees of freedom (df) is (8 categories - 1) = 7.")
    critical_value = scipy.stats.chi2.ppf(0.95, df=7)
    print(f"The critical value for a significance level of 0.05 and df=7 is {critical_value:.2f}.")
    print(f"Since {results[best_choice]:.2f} > {critical_value:.2f}, the null hypothesis is rejected.")
    
    # Output the final answer in the specified format
    print(f"\n<<<E>>>")

solve_chi_square()
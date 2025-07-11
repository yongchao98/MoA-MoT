import scipy.stats

def calculate_chi_square():
    """
    Calculates the chi-square statistic for different observed offspring phenotypes
    from a trihybrid cross, and determines which is most likely to reject the
    null hypothesis of independent assortment.
    """
    # Phenotypes in order:
    # T_R_Y_, T_R_yy, T_rrY_, T_rryy, ttR_Y_, ttR_yy, ttrrY_, ttrryy
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]
    
    # Expected ratio for a trihybrid cross with independent assortment
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    
    # Observed counts from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    results = {}
    
    print("--- Chi-Square Test for Independent Assortment ---")
    print(f"Null Hypothesis: Genes for height, seed shape, and seed color assort independently.")
    print(f"Expected phenotypic ratio: 27:9:9:3:9:3:3:1")
    # Degrees of freedom = 8 categories - 1 = 7
    # Significance level (alpha) = 0.05
    # Critical value for df=7, alpha=0.05 is approx 14.07
    print(f"Degrees of freedom: {len(phenotypes) - 1}")
    print(f"Significance level (alpha): 0.05")
    print(f"Critical Chi-Square Value: {scipy.stats.chi2.ppf(0.95, 7):.2f}\n")


    for choice, observed in observed_data.items():
        total_observed = sum(observed)
        if total_observed == 0:
            results[choice] = 0
            continue
        
        expected_counts = [(ratio/64) * total_observed for ratio in expected_ratio]
        
        chi_square_value = 0
        equation_terms = []
        for i in range(len(observed)):
            o = observed[i]
            e = expected_counts[i]
            if e == 0:
                # This case shouldn't happen here as all expected counts are > 0
                term = 0
            else:
                term = (o - e)**2 / e
            chi_square_value += term
            equation_terms.append(f"({o} - {e:.2f})^2 / {e:.2f}")

        results[choice] = chi_square_value
        
        print(f"--- Analysis for Choice {choice} ---")
        print(f"Observed counts: {observed}")
        print(f"Total Offspring: {total_observed}")
        print(f"Expected counts: {[round(e, 2) for e in expected_counts]}")
        print(f"Chi-Square Calculation: {' + '.join(equation_terms)} = {chi_square_value:.2f}")
        print("-" * (25 + len(choice)))
        print()
        
    # Find the choice with the maximum chi-square value
    max_choice = max(results, key=results.get)
    max_value = results[max_choice]
    
    print("\n--- Conclusion ---")
    print("The chi-square statistic measures the deviation from the expected ratio.")
    print("A higher value indicates a greater deviation and a higher likelihood of rejecting the null hypothesis.")
    for choice, value in results.items():
        print(f"Choice {choice}: Chi-Square = {value:.2f}")
    
    print(f"\nChoice {max_choice} has the highest chi-square value ({max_value:.2f}).")
    print("This extreme deviation from the expected 27:9:9:3:9:3:3:1 ratio makes it the most likely to lead to the rejection of the hypothesis of independent assortment.")


calculate_chi_square()
import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the Chi-square value for each answer choice to determine which
    is most likely to lead to the rejection of the null hypothesis of
    independent assortment in a trihybrid cross.
    """
    # Phenotype order:
    # 1. Tall, round, yellow (T_ R_ Y_)
    # 2. Tall, round, green (T_ R_ yy)
    # 3. Tall, wrinkled, yellow (T_ rr Y_)
    # 4. Tall, wrinkled, green (T_ rr yy)
    # 5. dwarf, round, yellow (tt R_ Y_)
    # 6. dwarf, round, green (tt R_ yy)
    # 7. dwarf, wrinkled, yellow (tt rr Y_)
    # 8. dwarf, wrinkled, green (tt rr yy)
    
    # Expected ratio for a TtRrYy x TtRrYy cross with independent assortment
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Observed counts for each answer choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    results = {}
    
    print("Calculating Chi-Square value for each option...")
    
    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        expected_counts = [(count / total_ratio_parts) * total_observed for count in expected_ratio]
        
        # Using scipy.stats.chisquare for accuracy
        # It calculates chi-square = sum((observed-expected)^2 / expected)
        chi2_statistic, p_value = scipy.stats.chisquare(f_obs=observed_counts, f_exp=expected_counts)
        
        results[option] = chi2_statistic
        print(f"Option {option}: Total Offspring = {total_observed}, Chi-Square = {chi2_statistic:.2f}")

    # Find the option with the highest Chi-Square value
    most_likely_rejection = max(results, key=results.get)
    
    print(f"\nThe option most likely to lead to rejection of the null hypothesis is '{most_likely_rejection}' because it has the highest Chi-Square value.")
    
    # Show the detailed calculation for the winning option
    print(f"\nDetailed calculation for Option {most_likely_rejection}:")
    
    winning_observed = observed_data[most_likely_rejection]
    total_observed = sum(winning_observed)
    winning_expected = [(count / total_ratio_parts) * total_observed for count in expected_ratio]
    
    calculation_string = "χ² = "
    terms = []
    for obs, exp in zip(winning_observed, winning_expected):
        terms.append(f"({obs} - {exp:.2f})² / {exp:.2f}")
    
    calculation_string += " + ".join(terms)
    
    print(calculation_string)
    print(f"Total Chi-Square Value = {results[most_likely_rejection]:.2f}")


solve_chi_square_problem()
<<<E>>>
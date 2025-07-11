import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square statistic for different observed outcomes of a trihybrid cross
    to find which is most likely to reject the null hypothesis of independent assortment.
    """
    phenotypes = [
        "tall, round, yellow", "tall, round, green", "tall, wrinkled, yellow",
        "tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]
    
    # Expected ratio for a Tt Rr Yy x Tt Rr Yy cross
    expected_ratios = [27, 9, 9, 3, 9, 3, 3, 1]
    
    # Observed counts from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    results = {}

    for option, observed_counts in observed_data.items():
        total_offspring = sum(observed_counts)
        
        # Calculate expected counts
        expected_counts = [(ratio / 64) * total_offspring for ratio in expected_ratios]
        
        # Calculate chi-square value
        chi_square_value = sum(
            ((o - e)**2) / e for o, e in zip(observed_counts, expected_counts)
        )
        results[option] = {
            "chi_square": chi_square_value,
            "observed": observed_counts,
            "expected": expected_counts
        }

    # Print the chi-square values for all options
    print("Chi-Square Test Results:")
    for option, data in results.items():
        print(f"Option {option}: χ² = {data['chi_square']:.2f}")

    # Find the option with the maximum chi-square value
    max_option = max(results, key=lambda k: results[k]['chi_square'])
    max_data = results[max_option]

    print(f"\nOption {max_option} has the highest chi-square value and is most likely to lead to rejection of the null hypothesis.")
    
    # Print the detailed equation for the max option
    print(f"\nDetailed calculation for Option {max_option}:")
    equation_parts = []
    for o, e in zip(max_data['observed'], max_data['expected']):
        equation_parts.append(f"({o} - {e:.2f})² / {e:.2f}")
    
    full_equation = " + ".join(equation_parts)
    print(f"χ² = {full_equation}")
    print(f"χ² = {max_data['chi_square']:.2f}")
    
    # Get the critical value for comparison
    df = len(phenotypes) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)
    print(f"\nThe degrees of freedom (df) is {df}. At a significance level of {alpha}, the critical value is {critical_value:.2f}.")
    print(f"Since {max_data['chi_square']:.2f} > {critical_value:.2f}, the null hypothesis would be rejected for Option {max_option}.")

solve_chi_square()
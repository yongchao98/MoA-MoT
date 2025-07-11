import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square value for each answer choice to determine which
    is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    # Expected phenotypic ratio for a trihybrid cross (27:9:9:3:9:3:3:1)
    ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio = sum(ratio_parts)

    # Observed counts for each answer choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    # Degrees of freedom = number of categories - 1
    df = len(ratio_parts) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print(f"The null hypothesis assumes a {':'.join(map(str, ratio_parts))} phenotypic ratio.")
    print(f"With {df} degrees of freedom and a significance level of {alpha}, the critical Chi-Square value is {critical_value:.2f}.")
    print("A calculated value greater than this will lead to rejection of the null hypothesis.")
    print("-" * 60)

    results = {}

    # Calculate Chi-Square for each choice
    for choice, observed in observed_data.items():
        total_observed = sum(observed)
        expected = [(part / total_ratio) * total_observed for part in ratio_parts]
        chi_square_value = sum([((o - e)**2) / e for o, e in zip(observed, expected)])
        results[choice] = chi_square_value
        print(f"Choice {choice}: Total Offspring = {total_observed}, Calculated Chi-Square = {chi_square_value:.2f}")

    # Find the choice with the maximum chi-square value
    max_choice = max(results, key=results.get)
    
    print("-" * 60)
    print(f"The combination most likely to lead to rejection is the one with the highest Chi-Square value.")
    print(f"The highest value is from Choice {max_choice} (χ² = {results[max_choice]:.2f}).")
    print("\nDetailed calculation for the final answer (Choice E):")
    
    # Print the detailed equation for the winning choice
    observed = observed_data[max_choice]
    total_observed = sum(observed)
    expected = [(part / total_ratio) * total_observed for part in ratio_parts]
    
    print("χ² = ", end="")
    final_terms = []
    for i in range(len(observed)):
        O = observed[i]
        E = expected[i]
        term_str = f"(({O} - {E:.2f})^2) / {E:.2f}"
        if i > 0:
            print(" + ", end="")
        print(term_str, end="")
        final_terms.append(((O - E)**2) / E)
    print()
    
    print("   = ", end="")
    for i, term in enumerate(final_terms):
      if i > 0:
        print(" + ", end="")
      print(f"{term:.2f}", end="")
    print()

    print(f"   = {sum(final_terms):.2f}")


solve_chi_square_problem()
<<<E>>>
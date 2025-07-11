def solve_chi_square_problem():
    """
    Calculates the chi-square value for each set of offspring phenotypes to determine
    which is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    
    # Phenotypes are ordered as they appear in the question.
    # The corresponding expected ratio parts for a trihybrid cross (27:9:9:3:9:3:3:1).
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]

    # Observed counts from each answer choice
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}
    
    print("--- Chi-Square Calculation Summary ---\n")

    # Store components for the best option for a detailed final printout
    best_option_details = {
        'option': None,
        'chi_square_value': -1,
        'terms': [],
        'equation': ''
    }

    # Loop through each option to calculate its chi-square value
    for option, observed in observed_counts.items():
        total_observed = sum(observed)
        total_ratio = sum(expected_ratio_parts)
        
        expected = [(part / total_ratio) * total_observed for part in expected_ratio_parts]
        
        chi_square_terms = []
        for i in range(len(observed)):
            O = observed[i]
            E = expected[i]
            # The chi-square test is not reliable if expected counts are zero.
            # In this problem, all expected counts are > 0.
            term = ((O - E)**2) / E
            chi_square_terms.append(term)

        chi_square_value = sum(chi_square_terms)
        results[option] = chi_square_value
        print(f"Option {option}: Chi-Square = {chi_square_value:.2f}")

        # Check if this is the best option so far
        if chi_square_value > best_option_details['chi_square_value']:
            best_option_details['option'] = option
            best_option_details['chi_square_value'] = chi_square_value
            best_option_details['terms'] = chi_square_terms
    
    # --- Analysis and Final Answer ---
    best_option = best_option_details['option']
    best_value = best_option_details['chi_square_value']
    
    print("\n--- Analysis ---")
    print("A higher Chi-Square value indicates a greater deviation from the expected ratio.")
    print(f"Option '{best_option}' has the highest Chi-Square value ({best_value:.2f}), meaning its observed results are the most different from what was expected under independent assortment.")
    print(f"Therefore, option '{best_option}' would most likely lead to the rejection of the null hypothesis.\n")

    print(f"--- Detailed Final Equation for Option {best_option} ---")
    print("The final Chi-Square value is the sum of each calculated term [(Observed - Expected)^2 / Expected]:\n")
    
    equation_terms_str = [f"{term:.2f}" for term in best_option_details['terms']]
    final_equation_str = " + ".join(equation_terms_str)
    
    print(f"Chi-Square = {final_equation_str}")
    print(f"Total Chi-Square = {best_value:.2f}")

solve_chi_square_problem()
<<<E>>>
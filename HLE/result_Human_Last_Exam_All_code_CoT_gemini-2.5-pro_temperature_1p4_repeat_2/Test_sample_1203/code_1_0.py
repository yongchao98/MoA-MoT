import sys

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for five different sets of genetic cross data
    to determine which is most likely to reject the hypothesis of independent assortment.
    """
    # The expected phenotypic ratio for a Tt Rr Yy x Tt Rr Yy cross is 27:9:9:3:9:3:3:1.
    # The phenotypes are ordered as:
    # 1. tall, round, yellow
    # 2. tall, round, green
    # 3. tall, wrinkled, yellow
    # 4. tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Observed counts from the answer choices
    options = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}

    for name, observed_counts in options.items():
        total_observed = sum(observed_counts)
        
        # Prevent division by zero if there are no offspring
        if total_observed == 0:
            print(f"\n--- Option {name}: No offspring observed. ---")
            results[name] = 0
            continue

        expected_counts = [(r / total_ratio_parts) * total_observed for r in expected_ratio]
        
        chi_square_value = 0
        equation_terms = []

        for i in range(len(observed_counts)):
            O = observed_counts[i]
            E = expected_counts[i]
            
            # The Chi-square test is not reliable if expected counts are very low (e.g., <5),
            # but we will calculate it anyway for this problem.
            if E == 0:
                 # This case won't happen here, but it's a good practice to handle it.
                term_value = float('inf') if O > 0 else 0
            else:
                term_value = ((O - E)**2) / E

            chi_square_value += term_value
            # For the equation, O is the observed count, E is the expected count.
            equation_terms.append(f"(({O} - {E:.2f})^2 / {E:.2f})")

        print(f"\n--- Analysis for Option {name} ---")
        equation_string = " + ".join(equation_terms)
        print(f"Chi-Square Equation = {equation_string}")
        print(f"Result: Chi-Square (χ²) for Option {name} = {chi_square_value:.2f}")
        results[name] = chi_square_value
        
    # Find the option with the maximum chi-square value
    # The `max` function with a `key` argument finds the item in a dictionary with the maximum value.
    most_likely_rejection = max(results, key=results.get)
    
    print("\n--- Conclusion ---")
    print("A higher chi-square value indicates a greater deviation from the expected ratio based on the hypothesis of independent assortment.")
    print(f"Comparing the results, Option {most_likely_rejection} has the highest chi-square value of {results[most_likely_rejection]:.2f}.")
    print("Therefore, this combination of offspring phenotypes would most likely lead to the rejection of the hypothesis.")

# Run the analysis
solve_chi_square_problem()
<<<E>>>
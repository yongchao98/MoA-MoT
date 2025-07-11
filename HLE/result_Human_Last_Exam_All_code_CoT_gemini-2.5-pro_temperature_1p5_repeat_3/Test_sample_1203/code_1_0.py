import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each of the given answer choices
    to determine which would most likely lead to the rejection of the null
    hypothesis of independent assortment.
    """
    # Phenotype categories in the standard 27:9:9:9:3:3:3:1 order
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Dwarf, round, yellow", "Tall, wrinkled, green", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]
    
    # Expected Mendelian ratios for a trihybrid cross
    ratios = [27, 9, 9, 9, 3, 3, 3, 1]
    total_ratio_parts = sum(ratios)

    # Observed counts from the answer choices.
    # Note: The order of phenotypes in the question differs from the standard ratio order.
    # We re-order the data here to match the 27:9:9:9:3:3:3:1 ratio order.
    # Question Order: TRY, TRg, TwY, Twg, dRY, dRg, dwY, dwg
    # Standard Ratio Order: TRY, TRg, TwY, dRY, Twg, dRg, dwY, dwg
    data = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],  # Swapping dRY and Twg to match ratio order
        "B": [180, 0, 0, 60, 0, 0, 0, 60],    # Swapping dRY and Twg
        "C": [144, 45, 45, 52, 16, 16, 16, 16],   # Swapping dRY and Twg
        "D": [150, 60, 50, 30, 40, 40, 30, 50],   # Swapping dRY and Twg
        "E": [0, 180, 0, 0, 0, 180, 0, 0]      # Swapping dRY and Twg
    }
    # Re-ordering the lists to match the ratio order: [TRY, TRg, TwY, dRY, Twg, dRg, dwY, dwg]
    # Standard order is: idx 0, 1, 2, 4, 3, 5, 6, 7 from question
    # This corresponds to ratio: 27, 9, 9, 9, 3, 3, 3, 1
    # We must match observed with correct ratio:
    # TRY (27), TRg (9), TwY (9), dRY(9), Twg(3), dRg(3), dwY(3), dwg(1)
    
    reordered_data = {}
    for key, values in data.items():
        reordered_data[key] = [
            values[0],  # TRY
            values[1],  # TRg
            values[2],  # TwY
            values[4],  # dRY
            values[3],  # Twg
            values[5],  # dRg
            values[6],  # dwY
            values[7]   # dwg
        ]
    
    results = {}

    for option, observed_counts in reordered_data.items():
        print(f"--- Analyzing Option {option} ---")
        
        total_observed = sum(observed_counts)
        
        # Calculate expected counts for each phenotype
        expected_counts = [(r / total_ratio_parts) * total_observed for r in ratios]
        
        chi_square_total = 0
        equation_parts = []
        
        print(f"Total Offspring = {total_observed}")
        print("Calculation: Sum of (Observed - Expected)^2 / Expected")

        for i in range(len(observed_counts)):
            o = observed_counts[i]
            e = expected_counts[i]
            # To avoid division by zero error, though not an issue here
            if e == 0:
                component = 0
            else:
                component = (o - e)**2 / e
            chi_square_total += component
            equation_parts.append(f"({o} - {e:.2f})^2 / {e:.2f}")

        # The problem asks to output each number in the final equation.
        # We will show the final calculated component for each phenotype.
        print("Chi-Square Equation Terms:")
        final_terms_str = []
        final_components = [(o - e)**2 / e for o, e in zip(observed_counts, expected_counts)]
        for term in final_components:
            final_terms_str.append(f"{term:.2f}")
        print(" + ".join(final_terms_str))

        print(f"Total Chi-Square Value (χ²) = {chi_square_total:.2f}")
        results[option] = chi_square_total
        print("-" * 25 + "\n")

    # Find and announce the winner
    most_likely_rejection = max(results, key=results.get)
    
    print("--- Conclusion ---")
    print("A higher Chi-Square value indicates a greater deviation from the expected ratio, making it more likely that the null hypothesis (independent assortment) is rejected.")
    print("Chi-Square values for each option:")
    for option, val in sorted(results.items()):
        print(f"Option {option}: {val:.2f}")
    
    print(f"\nOption {most_likely_rejection} has the highest Chi-Square value ({results[most_likely_rejection]:.2f}) and is therefore the most likely to lead to the rejection of the hypothesis of independent assortment.")
    
solve_chi_square_problem()
print("<<<E>>>")
import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_chi_square():
    """
    Calculates the chi-square statistic for each answer choice to determine which
    is most likely to lead to the rejection of the null hypothesis of independent assortment.
    """
    # Phenotype labels for clarity in the final output
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    # Observed counts for each option from the problem description
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected Mendelian ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    results = {}

    # Loop through each option to calculate its chi-square value
    for option, observed in observed_counts.items():
        total_observed = sum(observed)
        chi_square_sum = 0
        
        # Calculate expected values and the chi-square statistic
        for i in range(len(observed)):
            obs = observed[i]
            ratio_part = expected_ratio[i]
            expected = (ratio_part / total_ratio_parts) * total_observed
            
            # The chi-square contribution from this category
            if expected == 0:
                # This case implies the hypothesis predicts zero, which is not the case here.
                # If observed is also 0, contribution is 0. If not, it's infinite.
                term_val = 0 if obs == 0 else float('inf')
            else:
                term_val = ((obs - expected)**2) / expected
            
            chi_square_sum += term_val
        
        results[option] = chi_square_sum

    print("--- Chi-Square Calculation Results ---")
    for option, chi_square in results.items():
        print(f"Option {option}: χ² = {chi_square:.2f}")
    
    # Find the option with the highest chi-square value
    most_likely_rejection_option = max(results, key=results.get)
    
    print(f"\nOption {most_likely_rejection_option} has the highest chi-square value, indicating the largest deviation from the expected 27:9:9:3:9:3:3:1 ratio. This option would most likely lead to the rejection of the hypothesis of independent assortment.")

    # Display the detailed calculation for the winning option
    print(f"\n--- Detailed Calculation for Option {most_likely_rejection_option} ---")
    
    winning_observed = observed_counts[most_likely_rejection_option]
    total_observed_winning = sum(winning_observed)
    equation_parts = []
    
    print(f"Total Observed Offspring: {total_observed_winning}")
    print("χ² = Σ [ (Observed - Expected)² / Expected ]\n")

    for i in range(len(winning_observed)):
        obs = winning_observed[i]
        exp = (expected_ratio[i] / total_ratio_parts) * total_observed_winning
        equation_parts.append(f"(({obs} - {exp:.3f})^2 / {exp:.3f})")

    # Print the full equation with all its terms
    equation_string = " + \n       ".join(equation_parts)
    print(f"χ² = {equation_string}")
    
    final_chi_square_value = results[most_likely_rejection_option]
    print(f"\nχ² = {final_chi_square_value:.2f}")
    
    # Print the final answer in the required format
    print(f"\n<<<{most_likely_rejection_option}>>>")

solve_chi_square()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
import sys
# Redirect stdout to a variable to capture the final answer without printing it
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()


def calculate_and_find_max_chi_square():
    """
    Calculates the chi-square statistic for all choices and identifies the one
    with the maximum value, indicating the highest likelihood of rejecting the
    null hypothesis of independent assortment.
    """
    # Phenotypic ratio parts for a trihybrid cross (27:9:9:3:9:3:3:1)
    ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(ratio_parts)

    # Observed counts for each answer choice
    choices = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}

    # Calculate chi-square for each choice
    for label, observed_counts in choices.items():
        total_observed = sum(observed_counts)
        expected_counts = [(part / total_ratio_parts) * total_observed for part in ratio_parts]
        
        chi_square_value = 0
        for i in range(len(observed_counts)):
            O = observed_counts[i]
            E = expected_counts[i]
            # Avoid division by zero, though not expected here
            if E == 0:
                term = float('inf') if O > 0 else 0
            else:
                term = ((O - E)**2) / E
            chi_square_value += term
        results[label] = chi_square_value

    # Find the choice with the maximum chi-square value
    best_choice_label = max(results, key=results.get)
    max_chi_square_value = results[best_choice_label]

    # Now, print the detailed calculation for the winning choice
    print("A chi-square test compares observed offspring counts to those expected under a specific genetic hypothesis.")
    print("In this case, the hypothesis is independent assortment, which predicts a 27:9:9:3:9:3:3:1 phenotypic ratio.")
    print("A higher chi-square value indicates a greater deviation from the expected ratio and a higher likelihood of rejecting the hypothesis.")
    print("\nThe script calculates the chi-square value for all options and identifies the highest one.\n")
    print("--- Chi-Square Calculation for the Most Likely Rejected Hypothesis ---")
    
    # Recalculate values for the best choice for printing
    observed = choices[best_choice_label]
    total_observed = sum(observed)
    expected = [(part / total_ratio_parts) * total_observed for part in ratio_parts]
    
    print(f"Choice: {best_choice_label}")
    print(f"Observed Counts: {observed}")
    print(f"Total Observed: {total_observed}")
    print(f"Expected Counts: {[f'{e:.2f}' for e in expected]}")
    
    # Construct and print the final equation string
    equation_terms = []
    for i in range(len(observed)):
        term_str = f"({observed[i]} - {expected[i]:.2f})^2/{expected[i]:.2f}"
        equation_terms.append(term_str)
        
    print("\nFinal Equation:")
    print(f"Chi-Square = {' + '.join(equation_terms)}")

    print("\nResult:")
    print(f"Chi-Square for Choice {best_choice_label} = {max_chi_square_value:.2f}\n")
    print(f"The critical value for 7 degrees of freedom at a significance level of 0.05 is approximately 14.07.")
    print(f"Choice {best_choice_label} has the highest chi-square value, making it the most likely combination to lead to a rejection of the hypothesis of independent assortment.")
    
    # Print the final answer to be captured
    print(f"<<<{best_choice_label}>>>", file=sys.stdout)


calculate_and_find_max_chi_square()

# Restore original stdout and get captured output
sys.stdout = old_stdout
output = mystdout.getvalue()

# Find the final answer in the captured output
final_answer = output[output.find("<<<")+3:output.find(">>>")]
output_without_answer = output.replace(f'<<<{final_answer}>>>', '')

# Print the captured output without the answer part
print(output_without_answer.strip())
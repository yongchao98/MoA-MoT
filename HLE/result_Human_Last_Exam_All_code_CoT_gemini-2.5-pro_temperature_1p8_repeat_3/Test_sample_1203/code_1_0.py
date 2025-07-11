import sys
# Redirect stdout to a string buffer to capture all prints
# and then print the captured string at the end. This is a
# workaround for platforms that might truncate long outputs.
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each answer choice to determine which
    is most likely to lead to rejection of the null hypothesis of independent assortment.
    """
    # Phenotype names for clarity
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]
    # Expected Mendelian ratio for a trihybrid cross
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Observed counts from the answer choices
    observed_data = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}

    print("--- Chi-Square Test for Independent Assortment ---")
    print("Null Hypothesis (H₀): The genes assort independently.")
    print("Expected Ratio: 27:9:9:3:9:3:3:1\n")
    
    # Calculate chi-square for each choice
    for choice, observed_counts in observed_data.items():
        total_offspring = sum(observed_counts)
        
        # Calculate expected counts for each phenotype
        expected_counts = [(count / total_ratio_parts) * total_offspring for count in expected_ratio]
        
        # Calculate the chi-square value
        # Note: If an expected value is 0 (which doesn't happen here), it would cause a division by zero.
        chi_square_value = sum([((o - e)**2) / e for o, e in zip(observed_counts, expected_counts)])
        
        results[choice] = chi_square_value
        print(f"Choice {choice}: Total Offspring = {total_offspring}, Calculated Χ² = {chi_square_value:.2f}")

    # Determine the choice with the highest chi-square value
    max_choice = max(results, key=results.get)
    max_chi_square = results[max_choice]

    print("\n--- Conclusion ---")
    print("The combination of phenotypes that would most likely lead to rejection of the hypothesis is the one with the largest chi-square (Χ²) value, as it indicates the greatest deviation from the expected results.")
    print(f"\nThe highest calculated chi-square value belongs to Choice {max_choice} (Χ² = {max_chi_square:.2f}).")
    
    # Display the final equation for the winning choice
    print(f"\nFinal calculation for Choice {max_choice}:")
    
    obs_final = observed_data[max_choice]
    total_final = sum(obs_final)
    exp_final = [(count / total_ratio_parts) * total_final for count in expected_ratio]
    
    equation_parts = []
    sum_parts = []
    for o, e in zip(obs_final, exp_final):
        term_value = ((o-e)**2)/e
        equation_parts.append(f"(({o} - {e:.2f})² / {e:.2f})")
        sum_parts.append(f"{term_value:.2f}")

    # Printing the equation with numbers
    equation_str = " + ".join(equation_parts)
    print(f"Χ² = {equation_str}")
    # Printing the equation with intermediate sums
    sum_str = " + ".join(sum_parts)
    print(f"Χ² = {sum_str}")
    # Printing final value
    print(f"Χ² = {max_chi_square:.2f}")

# Execute the function and print the captured output
solve_chi_square_problem()
print(mystdout.getvalue())
sys.stdout = old_stdout

<<<E>>>
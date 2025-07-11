import sys
import io

# Redirect stdout to capture prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for different observed outcomes of a
    trihybrid cross to determine which would most likely lead to rejecting
    the hypothesis of independent assortment.
    """
    # Phenotype order from the problem statement
    phenotypes = [
        "tall, round, yellow", "tall, round, green", "tall, wrinkled, yellow",
        "tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    # Expected ratio based on independent assortment (TtRrYy x TtRrYy)
    # The order matches the phenotypes list
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    ratio_total = sum(expected_ratio)  # Should be 64

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
        total_observed = sum(observed_counts)
        
        # Calculate expected counts
        expected_counts = [(ratio / ratio_total) * total_observed for ratio in expected_ratio]
        
        # Calculate chi-square value
        chi_square_value = 0
        for obs, exp in zip(observed_counts, expected_counts):
            if exp == 0:
                # This case implies an impossible outcome under H0, chi-square would be infinite
                # but our expected counts are never zero.
                # If observed is non-zero and expected is zero, chi-square is infinite.
                chi_square_value = float('inf')
                break
            chi_square_value += ((obs - exp)**2) / exp
            
        results[option] = {
            'chi_square': chi_square_value,
            'observed': observed_counts,
            'expected': expected_counts
        }

    # Find the option with the maximum chi-square value
    best_option = max(results, key=lambda k: results[k]['chi_square'])
    
    print("The chi-square (χ²) test is used to compare observed genetic outcomes with expected outcomes to assess if the difference is statistically significant.")
    print("The null hypothesis (H₀) assumes the genes assort independently, following a 27:9:9:3:9:3:3:1 phenotypic ratio.")
    print("A higher χ² value indicates a larger deviation from the expected ratio. We need to find the option with the highest χ² value.\n")

    # Print results for all options for clarity
    for option, data in results.items():
        print(f"Option {option}: χ² = {data['chi_square']:.2f}")

    # Display the detailed calculation for the winning option
    winning_result = results[best_option]
    obs_list = winning_result['observed']
    exp_list = winning_result['expected']
    
    print(f"\nOption {best_option} has the highest χ² value and is the most likely to lead to a rejection of the null hypothesis.")
    
    # Building the equation string with each number
    equation_parts = []
    for o, e in zip(obs_list, exp_list):
        equation_parts.append(f"(({o} - {e:.2f})² / {e:.2f})")
    
    equation_string = " + ".join(equation_parts)
    print("\nDetailed calculation for Option {}:".format(best_option))
    print(f"χ² = {equation_string}")
    
    term_values = [((o - e)**2) / e for o, e in zip(obs_list, exp_list)]
    term_string = " + ".join([f"{t:.2f}" for t in term_values])
    print(f"   = {term_string}")
    print(f"   = {winning_result['chi_square']:.2f}")

    # Conclusion
    df = len(phenotypes) - 1
    p_value = 0.05
    critical_value = 14.07 # for df=7, p=0.05
    
    print(f"\nFor this test, the degrees of freedom (df) is {df} (8 categories - 1).")
    print(f"At a significance level of {p_value}, the critical χ² value is {critical_value}.")
    print(f"Since the calculated χ² value of {winning_result['chi_square']:.2f} is much greater than the critical value of {critical_value}, we would strongly reject the null hypothesis of independent assortment.")
    print(f"Therefore, the combination in option {best_option} would most likely lead to rejection.")
    
    return best_option

# Execute the function and capture the output
final_answer = solve_chi_square_problem()

# Get the captured output
output_string = captured_output.getvalue()

# Restore stdout
sys.stdout = old_stdout

# Print the captured output to the user
print(output_string)

# Print the final answer in the required format
print(f"<<<{final_answer}>>>")
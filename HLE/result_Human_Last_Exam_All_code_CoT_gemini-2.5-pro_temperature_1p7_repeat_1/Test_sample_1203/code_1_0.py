import sys
# Redirect stdout to a string to capture print output for the final answer
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each answer choice to determine
    which is most likely to lead to rejection of the null hypothesis.
    """
    # Phenotypes are listed in the order given by the problem statement.
    # The corresponding expected ratio for a trihybrid cross is 27:9:9:3:9:3:3:1
    ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(ratio_parts)  # Should be 64

    # Observed counts from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Degrees of freedom = number of categories - 1
    df = len(ratio_parts) - 1
    # For a significance level of 0.05 and df=7, the critical value is ~14.07
    critical_value = 14.07

    results = {}
    
    print(f"Chi-Square Test for Independent Assortment (Trihybrid Cross)")
    print(f"Expected Ratio: 27:9:9:3:9:3:3:1")
    print(f"Significance Level α = 0.05, Degrees of Freedom = {df}, Critical Value ≈ {critical_value}")
    print("-" * 70)

    for choice, observed in observed_data.items():
        total_offspring = sum(observed)
        
        # Calculate expected counts for this choice's total
        expected = [(part / total_ratio_parts) * total_offspring for part in ratio_parts]
        
        chi_square_terms = []
        equation_terms = []
        
        for i in range(len(observed)):
            O = observed[i]
            E = expected[i]
            
            # Avoid division by zero, though not expected in this problem
            term = ((O - E)**2) / E if E != 0 else 0
            chi_square_terms.append(term)
            equation_terms.append(f"({O} - {E:.2f})²/{E:.2f}")

        chi_square_value = sum(chi_square_terms)
        results[choice] = chi_square_value

        # Print the detailed calculation for this choice
        print(f"Analysis for Choice {choice}:")
        print(f"Total Observed = {total_offspring}")
        
        # Assemble and print the full equation string as requested
        equation_str = " + ".join(equation_terms)
        print(f"χ² = {equation_str}")
        print(f"χ² = {chi_square_value:.2f}")
        print("-" * 70)

    # Find and print the final conclusion
    most_likely_rejection = max(results, key=results.get)
    
    print("\nSummary of Chi-Square Values:")
    for choice, value in sorted(results.items()):
        print(f"  Choice {choice}: χ² = {value:.2f}")
    
    print("\nConclusion:")
    print(f"Choice {most_likely_rejection} has the highest chi-square value ({results[most_likely_rejection]:.2f}).")
    print("A higher chi-square value indicates a greater deviation from the expected results, making it the most likely to lead to a rejection of the null hypothesis of independent assortment.")
    
# Run the analysis
solve_chi_square_problem()

# Restore stdout and print the captured output
sys.stdout = old_stdout
captured_output = mystdout.getvalue()
print(captured_output)
print("<<<E>>>")
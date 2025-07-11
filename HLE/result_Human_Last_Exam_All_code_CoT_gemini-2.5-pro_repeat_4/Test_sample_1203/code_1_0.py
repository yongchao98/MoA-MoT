import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which
    is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype order:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. Dwarf, round, yellow
    # 6. Dwarf, round, green
    # 7. Dwarf, wrinkled, yellow
    # 8. Dwarf, wrinkled, green
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])

    # Observed data from answer choices
    observed_data = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    results = {}

    # Calculate Chi-Square for each choice
    for choice, observed in observed_data.items():
        total_offspring = np.sum(observed)
        expected_counts = (expected_ratio / np.sum(expected_ratio)) * total_offspring
        
        # To avoid division by zero in the unlikely case E is zero, add a small epsilon
        # but in this problem, all expected counts are non-zero.
        chi_square_terms = (observed - expected_counts)**2 / expected_counts
        chi_square_value = np.sum(chi_square_terms)
        
        results[choice] = {
            'chi_square': chi_square_value,
            'observed': observed,
            'expected': expected_counts,
            'terms': chi_square_terms
        }

    # Find the choice with the maximum Chi-Square value
    best_choice = max(results, key=lambda k: results[k]['chi_square'])
    
    print(f"The chi-square test is used to compare observed data with the expected 27:9:9:3:9:3:3:1 ratio.")
    print(f"The calculation is: Chi-Square = sum of ((Observed - Expected)^2 / Expected) for each phenotype.")
    print("\n--- Chi-Square Calculations ---")
    for choice, data in sorted(results.items()):
        print(f"Choice {choice}: Chi-Square = {data['chi_square']:.2f}")

    print("\n--- Detailed Calculation for the Best Choice ---")
    final_data = results[best_choice]
    observed = final_data['observed']
    expected = final_data['expected']
    terms = final_data['terms']
    
    print(f"Choice {best_choice} has the highest chi-square value, indicating the greatest deviation from the expected ratio.")
    print(f"This makes it the most likely to be rejected.\n")
    print(f"Observed Counts (O): {list(observed)}")
    print(f"Expected Counts (E): {[round(x, 2) for x in expected]}")
    
    equation_parts = []
    for i in range(len(observed)):
        part = f"(({observed[i]} - {expected[i]:.2f})^2 / {expected[i]:.2f})"
        equation_parts.append(part)
        
    print("\nChi-Square Equation:")
    print(" + ".join(equation_parts))

    term_values = " + ".join([f"{term:.2f}" for term in terms])
    print(f"\n= {term_values}")
    
    print(f"\n= {final_data['chi_square']:.2f}")

    # The degrees of freedom (df) = 8 categories - 1 = 7.
    # The critical value for p=0.05 at df=7 is 14.07.
    # All choices except C might have a high chi-square, but E is the highest.
    print(f"\nWith a chi-square value of {final_data['chi_square']:.2f}, which is far greater than the critical value (14.07 for p=0.05, df=7), we would strongly reject the null hypothesis of independent assortment for choice {best_choice}.")
    print(f"\n<<<E>>>")

solve_chi_square()
import sys

def solve_chi_square():
    """
    Calculates the chi-square statistic for each of the five answer choices
    to determine which is most likely to lead to the rejection of the
    null hypothesis of independent assortment.
    """
    # Phenotypes in order:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. Dwarf, round, yellow
    # 6. Dwarf, round, green
    # 7. Dwarf, wrinkled, yellow
    # 8. Dwarf, wrinkled, green
    
    # Expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)
    
    # Observed counts from the answer choices
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    results = {}
    
    print("Calculating Chi-Square values for each option:")
    print("-" * 50)
    
    # Calculate Chi-Square for each option
    for option, observed in observed_counts.items():
        total_offspring = sum(observed)
        
        # Calculate expected counts based on the total offspring
        expected = [(total_offspring * part) / total_ratio_parts for part in expected_ratio_parts]
        
        # Calculate the chi-square value and build the equation string
        chi_square_value = 0
        equation_string = []
        for i in range(len(observed)):
            o = observed[i]
            e = expected[i]
            # Avoid division by zero, although not expected here if total_offspring > 0
            if e == 0:
                # If expected is 0 and observed is not, deviation is infinite.
                # If both are 0, term is 0.
                term = float('inf') if o > 0 else 0
            else:
                term = (o - e)**2 / e
            chi_square_value += term
            equation_string.append(f"({o:.2f} - {e:.2f})^2 / {e:.2f}")

        # Store the result
        results[option] = chi_square_value
        
        # Print the detailed calculation for the current option
        print(f"For option {option}:")
        print(f"Observed Counts: {observed}")
        print(f"Total Offspring: {total_offspring}")
        expected_formatted = [f"{e:.2f}" for e in expected]
        print(f"Expected Counts: {expected_formatted}")
        
        full_equation = " + ".join(equation_string)
        print(f"Calculation: {full_equation} = {chi_square_value:.2f}\n")

    # Find the option with the maximum chi-square value
    most_likely_rejection = max(results, key=results.get)
    
    print("-" * 50)
    print("Summary of Chi-Square Values:")
    for option, value in results.items():
        print(f"Option {option}: {value:.2f}")
    
    print("\nConclusion:")
    print(f"A higher chi-square value indicates a greater deviation from the expected ratio.")
    print(f"The highest chi-square value is {results[most_likely_rejection]:.2f} for option {most_likely_rejection}.")
    print(f"Therefore, option {most_likely_rejection} would most likely lead to rejection of the hypothesis of independent assortment.")

# Execute the function
solve_chi_square()
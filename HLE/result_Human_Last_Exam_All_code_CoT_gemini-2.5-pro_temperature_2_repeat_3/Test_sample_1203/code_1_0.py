import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square statistic for different observed offspring phenotypes
    to determine which is most likely to reject the hypothesis of independent assortment.
    """
    # Phenotype order:
    # 1. Tall, round, yellow (TRY)
    # 2. Tall, round, green (TRG)
    # 3. Tall, wrinkled, yellow (TWY)
    # 4. Tall, wrinkled, green (TWG)
    # 5. dwarf, round, yellow (dRY)
    # 6. dwarf, round, green (dRG)
    # 7. dwarf, wrinkled, yellow (dWY)
    # 8. dwarf, wrinkled, green (dWG)
    
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected phenotypic ratio from a TtRrYy x TtRrYy cross
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    ratio_total = np.sum(expected_ratio)

    results = {}

    print("Calculating Chi-Square (χ²) values for each option:")
    print("-" * 50)
    
    for option, observed_counts in observed_data.items():
        observed_counts = np.array(observed_counts)
        total_offspring = np.sum(observed_counts)
        
        # Calculate expected counts based on the total offspring
        expected_counts = (expected_ratio / ratio_total) * total_offspring
        
        # Calculate chi-square value, handling cases where expected count is zero
        chi_square_terms = np.divide((observed_counts - expected_counts)**2, expected_counts, 
                                     out=np.zeros_like(expected_counts, dtype=float), 
                                     where=expected_counts!=0)
        
        total_chi_square = np.sum(chi_square_terms)
        results[option] = total_chi_square
        print(f"Option {option}: Total Offspring = {total_offspring}, χ² = {total_chi_square:.2f}")

    # Find the option with the maximum chi-square value
    winner = max(results, key=results.get)
    max_chi_square = results[winner]

    print("-" * 50)
    print(f"\nConclusion:")
    print(f"Option {winner} has the highest chi-square value (χ² = {max_chi_square:.2f}).")
    print("A higher chi-square value indicates a greater deviation from the expected results under the null hypothesis (independent assortment).")
    print(f"Therefore, the data in Option {winner} would most likely lead to the rejection of the hypothesis.\n")

    # Display the full calculation for the winning option
    print(f"Detailed Chi-Square Calculation for Option {winner}:")
    
    observed_winner = np.array(observed_data[winner])
    total_offspring_winner = np.sum(observed_winner)
    expected_winner = (expected_ratio / ratio_total) * total_offspring_winner
    
    equation_parts = []
    for obs, exp in zip(observed_winner, expected_winner):
        part = f"({obs} - {exp:.2f})²/ {exp:.2f}"
        equation_parts.append(part)

    print("χ² = " + " + ".join(equation_parts))
    
    # Print the value of each term in the sum
    term_values = []
    for obs, exp in zip(observed_winner, expected_winner):
        # handle case where exp=0
        if exp == 0:
            term_value = 0 # This case shouldn't happen with our expected ratios
        else:
            term_value = ((obs - exp)**2) / exp
        term_values.append(f"{term_value:.2f}")
    
    print("χ² = " + " + ".join(term_values))
    print(f"χ² = {max_chi_square:.2f}")

solve_chi_square()
<<<E>>>
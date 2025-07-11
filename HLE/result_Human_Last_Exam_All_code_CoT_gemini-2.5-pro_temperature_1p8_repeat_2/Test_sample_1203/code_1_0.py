import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square statistic for different sets of observed genetic cross data
    to determine which is most likely to reject the null hypothesis of independent assortment.
    """
    # Phenotypes in the standard order for a Tt Rr Yy cross
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    # The expected Mendelian ratio for a trihybrid cross (27:9:9:3:9:3:3:1)
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    
    # The observed counts from the answer choices
    observed_data = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    results = {}
    
    print("Calculating chi-square values for each option...\n")

    for option, observed in observed_data.items():
        total_offspring = np.sum(observed)
        
        # Calculate expected counts based on the total offspring and the Mendelian ratio
        expected_counts = (expected_ratio / np.sum(expected_ratio)) * total_offspring
        
        # Calculate the chi-square statistic
        chi_square = np.sum((observed - expected_counts)**2 / expected_counts)
        
        results[option] = chi_square
        print(f"Option {option}: Total Offspring = {total_offspring}, Chi-Square (χ²) = {chi_square:.2f}")

    # Find the option with the maximum chi-square value
    max_chi_square_option = max(results, key=results.get)
    
    print("\n------------------------------------------------------\n")
    print(f"The option with the highest chi-square value is '{max_chi_square_option}'.")
    print("A higher chi-square value indicates a larger deviation from the expected ratio,")
    print("making it the most likely to lead to the rejection of the null hypothesis.\n")

    # Display the detailed calculation for the winning option
    print(f"Detailed calculation for Option {max_chi_square_option}:")
    
    winning_observed = observed_data[max_chi_square_option]
    winning_total = np.sum(winning_observed)
    winning_expected = (expected_ratio / np.sum(expected_ratio)) * winning_total
    
    equation_parts = []
    for i in range(len(winning_observed)):
        o = winning_observed[i]
        e = winning_expected[i]
        part = f"({o} - {e:.2f})^2 / {e:.2f}"
        equation_parts.append(part)
        
    final_equation = " + ".join(equation_parts)
    final_chi_square_value = results[max_chi_square_option]

    print(f"\nFinal Equation (χ² = Σ [ (Observed - Expected)² / Expected ] ):")
    print(f"χ² = {final_equation}")
    print(f"   = {final_chi_square_value:.2f}")

    print("\n<<<" + max_chi_square_option + ">>>")

solve_chi_square()
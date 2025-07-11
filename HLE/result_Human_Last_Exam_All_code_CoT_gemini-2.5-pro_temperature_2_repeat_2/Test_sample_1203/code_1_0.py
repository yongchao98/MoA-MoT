import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square statistic for different observed outcomes of a
    three-gene cross to find which is most likely to reject the null hypothesis
    of independent assortment.
    """
    
    # The expected phenotypic ratio for a trihybrid cross is 27:9:9:9:3:3:3:1.
    expected_ratio = np.array([27, 9, 9, 9, 3, 3, 3, 1])
    ratio_proportions = expected_ratio / np.sum(expected_ratio)

    # Observed counts for each answer choice from the problem.
    observed_data = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }
    
    chi_square_values = {}

    print("Analyzing Chi-Square Values for Each Option\n")
    print("Null Hypothesis: Genes assort independently (27:9:9:9:3:3:3:1 ratio).")
    print("A higher Chi-Square value means a greater deviation from the expected ratio.\n")


    for option, observed in observed_data.items():
        print(f"--- Calculating for Option {option} ---")
        
        # 1. Calculate total number of offspring.
        total_offspring = np.sum(observed)
        
        # 2. Calculate expected counts based on the total.
        expected_counts = total_offspring * ratio_proportions
        
        # 3. Calculate the chi-square statistic.
        chi_square_statistic = np.sum((observed - expected_counts)**2 / expected_counts)
        chi_square_values[option] = chi_square_statistic
        
        print(f"Total Offspring = {total_offspring}")
        
        # 4. Display the full equation for transparency, including each number.
        print("Chi-Square Equation: Σ [ (Observed - Expected)² / Expected ]")
        equation_parts = []
        for i in range(len(observed)):
            # Format numbers to two decimal places for readability.
            obs_val = observed[i]
            exp_val = expected_counts[i]
            part = f"(({obs_val} - {exp_val:.2f})² / {exp_val:.2f})"
            equation_parts.append(part)
        
        # Joining all parts of the sum for the final equation string
        equation = " + ".join(equation_parts)
        print(f"χ² = {equation}")
        print(f"Result for Option {option}: Chi-Square Value = {chi_square_statistic:.2f}\n")

    # Find the option with the maximum chi-square value.
    max_chi_option = max(chi_square_values, key=chi_square_values.get)
    max_chi_value = chi_square_values[max_chi_option]

    print("--- Conclusion ---")
    print("The option with the largest chi-square value deviates most from the expected Mendelian ratio")
    print("and is therefore the most likely to lead to the rejection of the null hypothesis.\n")
    print(f"The highest chi-square value is {max_chi_value:.2f}, corresponding to Option {max_chi_option}.")

solve_chi_square()
<<<E>>>
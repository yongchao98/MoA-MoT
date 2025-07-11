import numpy as np

def calculate_and_print_chi_square():
    """
    Calculates the chi-square statistic for each of the five answer choices
    to determine which one deviates the most from the expected Mendelian ratio.
    """
    
    # Phenotype order:
    # 1. Tall, round, yellow; 2. Tall, round, green; 3. Tall, wrinkled, yellow; 4. Tall, wrinkled, green;
    # 5. dwarf, round, yellow; 6. dwarf, round, green; 7. dwarf, wrinkled, yellow; 8. dwarf, wrinkled, green.
    phenotype_labels = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow", "Tall, wrinkled, green",
        "dwarf, round, yellow", "dwarf, round, green", "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    # Expected ratio for a trihybrid cross with independent assortment
    expected_ratio_parts = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    total_ratio_parts = np.sum(expected_ratio_parts)

    # Observed counts for each answer choice
    observed_data = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    results = {}

    print("Calculating Chi-Square values to find the data set most likely to be rejected.\n")

    for option, observed_counts in observed_data.items():
        total_offspring = np.sum(observed_counts)
        
        # Calculate expected counts based on the total number of offspring
        expected_counts = (expected_ratio_parts / total_ratio_parts) * total_offspring
        
        # Calculate the Chi-Square statistic terms
        chi_square_terms = (observed_counts - expected_counts)**2 / expected_counts
        
        chi_square_value = np.sum(chi_square_terms)
        results[option] = chi_square_value
        
        print(f"--- Option {option} ---")
        print(f"Total Offspring: {total_offspring}")
        print("\nChi-Square Calculation Terms [(Observed - Expected)^2 / Expected]:")
        
        # Print each term of the calculation
        for i in range(len(observed_counts)):
            o = observed_counts[i]
            e = expected_counts[i]
            term = chi_square_terms[i]
            # The following line shows each number in the equation for this term
            print(f"  {i+1:2}. ({o:3} - {e:6.2f})^2 / {e:6.2f} = {term:8.2f}   ({phenotype_labels[i]})")
        
        print(f"\nTotal Chi-Square (χ²) for Option {option}: {chi_square_value:.2f}\n" + "-"*20)

    # Find the option with the highest Chi-Square value
    max_chi_square_option = max(results, key=results.get)
    max_chi_square_value = results[max_chi_square_option]

    print("\n--- Conclusion ---")
    print("A higher Chi-Square value indicates a greater deviation from the expected ratio,")
    print("making it more likely to reject the null hypothesis of independent assortment.")
    print(f"\nOption {max_chi_square_option} has the highest Chi-Square value ({max_chi_square_value:.2f}).")
    print(f"Therefore, the data in Option {max_chi_square_option} would most likely lead to rejection of the hypothesis.")

# Run the calculation and print the results
calculate_and_print_chi_square()
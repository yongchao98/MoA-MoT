import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square statistic for each of the five options to determine
    which would most likely lead to the rejection of the hypothesis of independent assortment.
    """
    # Define the expected phenotypic ratio for a Tt Rr Yy x Tt Rr Yy cross
    expected_ratio = np.array([27, 9, 9, 9, 3, 3, 3, 1])
    ratio_sum = np.sum(expected_ratio) # Should be 64

    # Phenotype labels for clear output
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Dwarf, round, yellow", "Tall, wrinkled, green", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Observed counts from the answer choices
    observed_counts_options = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    results = {}

    print("The chi-square test determines the deviation of observed from expected outcomes.")
    print("The null hypothesis of independent assortment predicts a 27:9:9:9:3:3:3:1 ratio.")
    print("We will calculate the chi-square value ( Σ[(O-E)²/E] ) for each option.")
    print("A higher chi-square value indicates a greater deviation and is more likely to be rejected.\n")

    for option, observed in observed_counts_options.items():
        total_offspring = np.sum(observed)
        expected_counts = (expected_ratio / ratio_sum) * total_offspring

        # To avoid division by zero if an expected count was 0 (not the case here)
        chi_square_terms = np.divide(
            np.square(observed - expected_counts),
            expected_counts,
            out=np.zeros_like(expected_counts, dtype=float),
            where=expected_counts!=0
        )
        chi_square_value = np.sum(chi_square_terms)
        results[option] = chi_square_value

        print(f"--- Calculation for Option {option} ---")
        print(f"Observed (O): {observed.tolist()}")
        print(f"Total Offspring (N): {total_offspring}")
        print(f"Expected (E): {[round(e, 2) for e in expected_counts]}")
        
        # Build and print the equation string with all numbers
        equation_str = "Chi-Square = "
        terms_str_list = []
        for i in range(len(observed)):
            o = observed[i]
            e = expected_counts[i]
            term_str = f"({o}-{e:.2f})^2/{e:.2f}"
            terms_str_list.append(term_str)
        equation_str += " + ".join(terms_str_list)
        print(equation_str)

        values_str = "= "
        values_str_list = [f"{term:.2f}" for term in chi_square_terms]
        values_str += " + ".join(values_str_list)
        print(values_str)
        
        print(f"= {chi_square_value:.2f}\n")

    # Find the option with the highest chi-square value
    most_likely_rejection = max(results, key=results.get)
    
    # Set the significance level and degrees of freedom to compare
    df = len(expected_ratio) - 1
    critical_value = 14.07 # for p=0.05 and df=7
    
    print("--- Conclusion ---")
    print(f"The critical value for a chi-square test with a significance level of 0.05 and {df} degrees of freedom is {critical_value}.")
    print(f"Any chi-square value greater than {critical_value} would lead to the rejection of the null hypothesis.")
    print(f"\nOption {most_likely_rejection} has the highest chi-square value ({results[most_likely_rejection]:.2f}).")
    print("This indicates the largest deviation from the expected ratio of independent assortment.")
    print(f"Therefore, the combination of phenotypes in option {most_likely_rejection} would most likely lead to the rejection of the hypothesis.")

solve_chi_square()
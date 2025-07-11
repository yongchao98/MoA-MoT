import numpy as np

def calculate_chi_square():
    """
    Calculates the chi-square statistic for each of the five answer choices
    to determine which would most likely lead to rejecting the hypothesis of
    independent assortment.
    """

    # Phenotypes corresponding to the ratio parts
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]
    
    # Expected phenotypic ratio for a trihybrid cross (TtRrYy x TtRrYy)
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    total_ratio_parts = np.sum(expected_ratio)

    # Observed offspring counts for each answer choice
    observed_data = {
        "A": np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        "B": np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        "C": np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        "D": np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        "E": np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    print("Analyzing Chi-Square Values for Each Option\n" + "="*50)
    print(f"The null hypothesis of independent assortment predicts a phenotypic ratio of 27:9:9:3:9:3:3:1.")
    print(f"Degrees of freedom = (Number of phenotypes) - 1 = {len(phenotypes)} - 1 = 7.")
    print("The critical value at a significance level of 0.05 is 14.07.")
    print("A calculated Chi-Square value > 14.07 leads to rejection of the null hypothesis.")
    print("We are looking for the largest Chi-Square value.\n")


    results = {}

    for option, observed_counts in observed_data.items():
        total_offspring = np.sum(observed_counts)
        
        # Calculate expected counts for each phenotype
        expected_counts = (expected_ratio / total_ratio_parts) * total_offspring
        
        # Calculate chi-square value
        # Note: Add a small epsilon to the denominator to avoid division by zero
        # if an expected count were 0, though it's not the case here.
        chi_square_terms = (observed_counts - expected_counts)**2 / (expected_counts + 1e-9)
        chi_square_value = np.sum(chi_square_terms)
        results[option] = chi_square_value

        print(f"--- Option {option} ---")
        print(f"Observed counts: {observed_counts.tolist()}")
        print(f"Total Offspring: {total_offspring}")
        print(f"Expected counts: {[round(c, 2) for c in expected_counts]}")
        
        # Outputting the numbers in the final equation as a sum of terms
        equation_str = " + ".join([f"({o}-{e:.2f})^2/{e:.2f}" for o, e in zip(observed_counts, expected_counts)])
        terms_str = " + ".join([f"{term:.2f}" for term in chi_square_terms])
        
        print(f"\nEquation: Chi-Square = {equation_str}")
        print(f"Contribution of each term: {terms_str}")
        print(f"Total Chi-Square for Option {option}: {chi_square_value:.2f}\n")
        print("-"*50)

    # Find the option with the maximum chi-square value
    most_likely_rejection = max(results, key=results.get)
    max_chi_square = results[most_likely_rejection]

    print("\n--- Conclusion ---")
    print("Comparison of Chi-Square values:")
    for option, value in results.items():
        print(f"Option {option}: {value:.2f}")

    print(f"\nOption {most_likely_rejection} has the highest Chi-Square value ({max_chi_square:.2f}).")
    print("Therefore, this combination of offspring phenotypes would most likely lead to the rejection of the hypothesis of independent assortment.")

if __name__ == '__main__':
    calculate_chi_square()
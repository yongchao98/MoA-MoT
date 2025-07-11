def find_max_cutoff():
    """
    This script finds the maximum Aortomesenteric (AM) diameter cutoff for identifying
    a condition (like EVP enhancement, often associated with Nutcracker Syndrome)
    that meets specific sensitivity and specificity criteria.

    The analysis is based on data adapted from a medical study where the diagnostic
    test is defined as 'AM distance < cutoff'.
    """
    # Data from Gulleroglu et al. (2014) on AM distance for diagnosing Nutcracker Syndrome.
    # The list contains dictionaries, each with a cutoff threshold and its performance.
    data = [
        {'cutoff': 8, 'sensitivity': 95.8, 'specificity': 37.5},
        {'cutoff': 7, 'sensitivity': 83.3, 'specificity': 50.0},
        {'cutoff': 6, 'sensitivity': 79.2, 'specificity': 75.0},
        {'cutoff': 5, 'sensitivity': 75.0, 'specificity': 87.5},
        {'cutoff': 4, 'sensitivity': 50.0, 'specificity': 100.0},
        {'cutoff': 3, 'sensitivity': 16.7, 'specificity': 100.0}
    ]

    # Define the required sensitivity and specificity.
    sensitivity_requirement = 60
    specificity_requirement = 80

    qualifying_cutoffs = []

    print("Analyzing Aortomesenteric diameter cutoffs...")
    print(f"Criteria: Sensitivity > {sensitivity_requirement}% AND Specificity > {specificity_requirement}%")
    print("-" * 65)

    # Iterate through each data point to check if it meets the criteria.
    for record in data:
        cutoff = record['cutoff']
        sensitivity = record['sensitivity']
        specificity = record['specificity']

        # Check if the sensitivity and specificity are greater than the requirements.
        is_sensitive = sensitivity > sensitivity_requirement
        is_specific = specificity > specificity_requirement

        print(f"Evaluating test 'AM Diameter < {cutoff} mm':")
        # Print the comparison for sensitivity, including the numbers.
        print(f"  - Sensitivity Check: {sensitivity}% > {sensitivity_requirement}%  ->  {'Pass' if is_sensitive else 'Fail'}")
        # Print the comparison for specificity, including the numbers.
        print(f"  - Specificity Check: {specificity}% > {specificity_requirement}%  ->  {'Pass' if is_specific else 'Fail'}")

        # If both conditions are met, add the cutoff to the list of solutions.
        if is_sensitive and is_specific:
            qualifying_cutoffs.append(cutoff)
            print("  --> Result: This cutoff qualifies.\n")
        else:
            print("  --> Result: This cutoff does not qualify.\n")

    print("-" * 65)

    # Determine the final answer.
    if not qualifying_cutoffs:
        print("Conclusion: No AM diameter cutoff value in the dataset satisfies both conditions.")
    else:
        # Find the maximum value from the list of qualifying cutoffs.
        max_qualifying_cutoff = max(qualifying_cutoffs)
        print(f"The qualifying cutoff threshold(s) are: {qualifying_cutoffs} mm.")
        print(f"The maximum (at most) qualifying AM diameter threshold is {max_qualifying_cutoff} mm.")


if __name__ == "__main__":
    find_max_cutoff()
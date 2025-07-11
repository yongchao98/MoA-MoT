import sys

def calculate_simpsons_index():
    """
    This function demonstrates the calculation of Simpson's Diversity Index
    for a sample containing only one species, showing how a value of 0 is obtained.
    """
    # Step 1: Define the student's sample data.
    # The sample only contains individuals from one species.
    student_sample = {"Common Pipistrelle": 50}
    
    print("Scenario: A sample is collected containing 50 individuals of a single bat species.")
    print(f"Sample Data: {student_sample}\n")

    # Step 2: Calculate N, the total number of individuals.
    N = sum(student_sample.values())
    
    # Step 3: Calculate the sum of n(n-1) for each species.
    # Since there's only one species, this is straightforward.
    n = list(student_sample.values())[0]
    sum_n_n_minus_1 = n * (n - 1)
    
    # Step 4: Calculate the denominator N(N-1).
    N_N_minus_1 = N * (N - 1)

    print("---[Mathematical Calculation]---")
    print(f"The formula for Simpson's Diversity Index is D = 1 - [Σn(n-1) / N(N-1)]")
    
    # Check to avoid division by zero if N < 2
    if N < 2:
        print("\nCannot calculate Simpson's Index for a sample with fewer than 2 individuals.")
        return

    # Step 5: Calculate Simpson's Diversity Index (D).
    dominance_index = sum_n_n_minus_1 / N_N_minus_1
    D = 1 - dominance_index

    # Step 6: Print the final equation with all numbers.
    print("\nSubstituting the values into the formula:")
    print(f"Total individuals (N) = {N}")
    print(f"Individuals of the single species (n) = {n}")
    print(f"Value of n(n-1) = {n}({n-1}) = {sum_n_n_minus_1}")
    print(f"Value of N(N-1) = {N}({N-1}) = {N_N_minus_1}")
    
    # We must explicitly show each number in the final equation.
    print("\nFinal Equation:")
    print(f"D = 1 - ( {sum_n_n_minus_1} / {N_N_minus_1} )")
    print(f"D = 1 - {dominance_index}")
    print(f"D = {D}\n")

    # Step 7: Provide the final interpretation.
    print("---[Conclusion]---")
    print("1. Mathematical Validity: The calculation shows that an index of 0 is a valid mathematical result when the sample contains only one species.")
    print("2. Ecological Validity: The result (D=0) suggests zero diversity. However, it is known that the island has multiple bat species. Therefore, the sample is not representative of the true ecological reality of the island.")
    print("\nThe index is mathematically valid, but not ecologically valid.")

# Run the function
calculate_simpsons_index()
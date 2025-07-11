def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # Calculate powers of q
    q2 = q**2
    q3 = q**3
    q4 = q**4

    # Calculate factors for the first term
    term1_factor1 = q2 - q + 1
    term1_factor2 = q2 + 1
    
    # Calculate the first term
    term1 = q4 * term1_factor1 * term1_factor2

    # Calculate factors for the second term
    term2_factor1 = q3 + 1
    
    # Calculate the second term
    term2 = q2 * term2_factor1
    
    # Calculate the total number of involutions
    total_involutions = term1 + term2

    # Print the step-by-step evaluation of the formula
    print("The formula for the number of involutions in PSU(4, q) is:")
    print("N = q^4 * (q^2 - q + 1) * (q^2 + 1) + q^2 * (q^3 + 1)")
    print("\nFor q = 997:")
    print(f"N = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1) + {q}^2 * ({q}^3 + 1)")
    print("\nEvaluating the powers of q:")
    print(f"q^2 = {q2}")
    print(f"q^3 = {q3}")
    print(f"q^4 = {q4}")
    print("\nSubstituting these values into the formula:")
    print(f"N = {q4} * ({q2} - {q} + 1) * ({q2} + 1) + {q2} * ({q3} + 1)")
    print("\nEvaluating the expressions in parentheses:")
    print(f"N = {q4} * ({term1_factor1}) * ({term1_factor2}) + {q2} * ({term2_factor1})")
    print("\nCalculating the two main terms:")
    print(f"First term  = {q4} * {term1_factor1} * {term1_factor2} = {term1}")
    print(f"Second term = {q2} * {term2_factor1} = {term2}")
    print("\nAdding the two terms to get the final result:")
    print(f"Total number of involutions = {term1} + {term2} = {total_involutions}")

solve()
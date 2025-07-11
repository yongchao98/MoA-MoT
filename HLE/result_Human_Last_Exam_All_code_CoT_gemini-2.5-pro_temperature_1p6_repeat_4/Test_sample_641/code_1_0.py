def solve_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions is given by the formula:
    # N = q^4 * (q + 1) * (q^3 + 1) * (q^2 + 1)
    
    # Calculate each term
    term1 = q**4
    term2 = q + 1
    term3 = q**3 + 1
    term4 = q**2 + 1

    # Calculate the final result
    result = term1 * term2 * term3 * term4

    # Print the equation and the final answer
    print("The number of involutions in PSU(4,997) is given by the equation:")
    print(f"N = q^4 * (q + 1) * (q^3 + 1) * (q^2 + 1) for q = 997")
    print("Substituting the value of q:")
    print(f"N = {term1} * {term2} * {term3} * {term4}")
    print("\nFinal Answer:")
    print(result)

solve_involutions()
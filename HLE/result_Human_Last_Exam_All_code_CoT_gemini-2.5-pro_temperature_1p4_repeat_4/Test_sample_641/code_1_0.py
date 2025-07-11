def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions is:
    # (1/2) * q^4 * (q^2 - q + 1) * (q^2 + 1) * (q^3 + 2q^2 + q + 1)
    
    # Calculate each term of the formula
    term1 = q**4
    term2 = q**2 - q + 1
    term3 = q**2 + 1
    term4 = q**3 + 2*q**2 + q + 1

    # The total number is the product of the terms divided by 2
    # We use integer division // as the result is guaranteed to be an integer.
    num_involutions = (term1 * term2 * term3 * term4) // 2
    
    # Print the equation with the computed values of each term
    print("The number of involutions is given by the equation:")
    print(f"({term1} * {term2} * {term3} * {term4}) / 2")
    print("\nWhich evaluates to:")
    print(num_involutions)

solve()
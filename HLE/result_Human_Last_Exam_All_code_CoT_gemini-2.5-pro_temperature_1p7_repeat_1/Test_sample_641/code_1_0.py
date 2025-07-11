def solve():
    """
    This function calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions is given by the formula:
    # Num = q^4 * (q^4 - q^3 + q^2 - q + 1)
    # We calculate the two parts of the product separately.

    term1 = q**4
    term2 = q**4 - q**3 + q**2 - q + 1
    
    result = term1 * term2
    
    # Print the final equation with the computed numbers.
    print(f"The number of involutions in PSU(4, 997) is given by the expression q^4 * (q^4 - q^3 + q^2 - q + 1) for q=997.")
    print(f"The final equation with computed values is:")
    print(f"{term1} * {term2} = {result}")

solve()
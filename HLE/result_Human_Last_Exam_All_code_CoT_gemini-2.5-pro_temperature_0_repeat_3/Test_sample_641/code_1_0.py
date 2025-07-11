def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions in PSU(4, q) for q odd
    # is (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2.

    # Calculate the terms of the formula
    q2 = q * q
    q4 = q2 * q2
    term1 = q2 - q + 1
    term2 = q2 + 1

    # Calculate the total number of involutions
    # Using integer division // as the result is guaranteed to be an integer.
    num_involutions = (q4 * term1 * term2) // 2

    # Print the equation with the numbers substituted
    print("The number of involutions is calculated using the formula:")
    print("(q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2")
    print("\nSubstituting q = 997:")
    print(f"({q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1)) / 2")
    print(f"= ({q4} * ({q2} - {q} + 1) * ({q2} + 1)) / 2")
    print(f"= ({q4} * {term1} * {term2}) / 2")
    print(f"= {q4 * term1 * term2} / 2")
    print(f"= {num_involutions}")

solve()
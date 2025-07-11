def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    This corresponds to the order of the general linear group GL(4, 2).
    """
    # The parameters are n=4 and q=2 for GL(n, q).
    n = 4
    q = 2

    # The formula for the order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.
    term1 = q**n - q**0  # (16 - 1)
    term2 = q**n - q**1  # (16 - 2)
    term3 = q**n - q**2  # (16 - 4)
    term4 = q**n - q**3  # (16 - 8)

    # Calculate the total order
    result = term1 * term2 * term3 * term4

    # Print the final equation with each number as requested
    print("The highest possible order for E is |GL(4, 2)|.")
    print("The calculation is:")
    print(f"{term1} * {term2} * {term3} * {term4} = {result}")

solve()
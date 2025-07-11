def solve():
    """
    Calculates the number of internal adjunctions from F_11^3 to itself
    in the described 2-category. This corresponds to the order of the
    general linear group GL(3, F_11).
    """
    q = 11
    n = 3

    # The order of GL(n, F_q) is the product of (q^n - q^i) for i from 0 to n-1.
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2

    # Calculate the final product
    result = term1 * term2 * term3

    # Print the equation with the calculated numbers, as requested.
    print("The number of adjunctions is calculated by the following equation:")
    print(f"{term1} * {term2} * {term3} = {result}")

solve()
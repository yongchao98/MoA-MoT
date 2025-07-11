def calculate_norm_for_Jn(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn,
    where n is an odd positive integer.

    The derivation shows that the 1-norm is given by the formula:
    ||T||_1 = 4**(n+1) - 2**(n+1)

    Args:
        n (int): An odd positive integer.

    Returns:
        int: The 1-norm of the correlation matrix T, or None if n is invalid.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print(f"Error: n must be an odd positive integer. Received n = {n}.")
        return None

    # The derived formula for the 1-norm
    power = n + 1
    term1 = 4**power
    term2 = 2**power
    result = term1 - term2

    # Output the result with the equation details
    print(f"For n = {n}:")
    print(f"The 1-norm is calculated using the formula: 4^(n+1) - 2^(n+1)")
    print(f"4^({n}+1) - 2^({n}+1) = 4^{power} - 2^{power} = {term1} - {term2} = {result}")
    print("-" * 30)

    return result

if __name__ == '__main__':
    # Demonstrate the calculation for a few odd values of n
    calculate_norm_for_Jn(1)
    calculate_norm_for_Jn(3)
    calculate_norm_for_Jn(5)
    
    # Example of invalid input
    calculate_norm_for_Jn(2)

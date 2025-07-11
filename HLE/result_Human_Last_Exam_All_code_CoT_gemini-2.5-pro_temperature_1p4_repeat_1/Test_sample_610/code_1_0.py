def solve_l_exact(n):
    """
    Calculates the exact value of l(n, b) based on the symbolic derivation.
    The result is 2n(n-1).

    Args:
    n (int): An integer >= 10.

    Returns:
    int: The calculated value of l(n, b).
    """
    if not isinstance(n, int) or n < 10:
        print("Error: n must be an integer greater than or equal to 10.")
        return None
    
    # Perform the calculation based on the derived formula
    n_minus_1 = n - 1
    result = 2 * n * n_minus_1
    
    # The final answer needs to output the full equation with numbers.
    print(f"2 * {n} * {n_minus_1} = {result}")
    return result

# Example usage with a value for n (e.g., n = 20, which satisfies n >= 10)
n_example = 20
solve_l_exact(n_example)
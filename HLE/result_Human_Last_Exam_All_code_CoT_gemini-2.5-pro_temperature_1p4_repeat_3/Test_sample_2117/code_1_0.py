def find_least_upper_bound():
    """
    This function calculates the least upper bound of the expression
    E_P * E_H * S_P * S_H.

    The problem is solved analytically, and the final expression for the product is
    a function of n: f(n) = 2 * (n + 2) / (n + 1).

    This function is decreasing for n >= 1. Its least upper bound is the value at n=1.
    This script calculates and prints this value.
    """
    # The least upper bound occurs at the smallest positive integer n, which is 1.
    n = 1
    
    # Calculate the components of the final equation
    term1 = 2
    term2_numerator = n + 2
    term2_denominator = n + 1

    # Calculate the final result
    result = term1 * term2_numerator / term2_denominator

    # As per the instructions, print each number in the final equation.
    print("The final expression for the product as a function of n is: 2 * (n + 2) / (n + 1)")
    print("This function is decreasing for positive n, so its least upper bound is the value for n = 1.")
    print(f"Calculation for n = {n}:")
    print(f"{term1} * ({n} + 2) / ({n} + 1) = {result}")

find_least_upper_bound()
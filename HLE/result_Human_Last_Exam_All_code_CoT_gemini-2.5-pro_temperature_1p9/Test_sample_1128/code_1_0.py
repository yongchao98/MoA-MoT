import math

def count_cool_strings(n):
    """
    Calculates the number of cool strings of maximal length for n symbols.
    The formula used is n! * (n+1)!.
    """
    if n < 1:
        print("Please enter a positive integer for n.")
        return

    n_factorial = math.factorial(n)
    n_plus_1_factorial = math.factorial(n + 1)
    
    result = n_factorial * n_plus_1_factorial
    
    # Printing the equation and the final result
    print(f"For n = {n}:")
    print(f"The number of cool strings is n! * (n+1)!")
    print(f"{n}! * {n+1}! = {n_factorial} * {n_plus_1_factorial} = {result}")

# Example usage with n=4, as described in the plan.
# The user can change this value to calculate for a different 'n'.
n_example = 4
count_cool_strings(n_example)
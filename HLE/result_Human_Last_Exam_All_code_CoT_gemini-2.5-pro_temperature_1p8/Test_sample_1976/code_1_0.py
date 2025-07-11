import math

def calculate_T_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for odd n.

    Args:
        n: An odd integer.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: n must be a positive odd integer.")
        return

    # According to the derived formula for odd n: ||T||_1 = (4^((n+1)/2) - 1) * 4^n
    
    # Calculate components of the formula
    base1 = 4
    exponent1 = (n + 1) // 2
    term1 = base1**exponent1
    
    factor1 = term1 - 1
    
    base2 = 4
    exponent2 = n
    factor2 = base2**exponent2
    
    # Calculate the final result
    result = factor1 * factor2
    
    # Print the result with the equation
    print(f"For n = {n}:")
    print(f"The 1-norm of the correlation matrix T is calculated by the formula: (4^((n+1)/2) - 1) * 4^n")
    print(f"Substituting n = {n}:")
    print(f"({base1}^(({n}+1)/2) - 1) * {base2}^{n}")
    print(f"= ({base1}^{exponent1} - 1) * {factor2}")
    print(f"= ({term1} - 1) * {factor2}")
    print(f"= {factor1} * {factor2}")
    print(f"= {result}")


# Example usage for n=3, a positive odd integer.
# You can change this value to any other positive odd integer.
n = 3
calculate_T_norm(n)

# Let's provide another example for n=5
# n = 5
# calculate_T_norm(n)
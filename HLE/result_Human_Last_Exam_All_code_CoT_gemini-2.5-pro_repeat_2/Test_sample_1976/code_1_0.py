import math

def calculate_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn with odd n.

    Args:
        n (int): An odd integer.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: Please provide a positive odd integer for n.")
        return

    # The formula for the 1-norm is 2**(n+1) * (2**(n+1) - 1)
    
    # Calculate the exponent
    exponent = n + 1
    
    # Calculate the first term, 2**(n+1)
    term1 = 2**exponent
    
    # Calculate the second term, (2**(n+1) - 1)
    term2 = term1 - 1
    
    # Calculate the final result
    result = term1 * term2
    
    # Print the equation with the calculated numbers
    print(f"For n = {n}:")
    print(f"The 1-norm of the correlation matrix is calculated by the formula: 2^(n+1) * (2^(n+1) - 1)")
    print(f"Substituting n = {n}, we get:")
    print(f"2^({n}+1) * (2^({n}+1) - 1) = 2^{exponent} * (2^{exponent} - 1)")
    print(f"= {term1} * ({term1} - 1)")
    print(f"= {term1} * {term2}")
    print(f"= {result}")

# --- User Input ---
# You can change the value of n here to any positive odd integer.
# For example, n = 1, 3, 5, etc.
n = 3

# Execute the calculation
calculate_norm(n)

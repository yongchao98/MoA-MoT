import math

def calculate_equilibria(N):
    """
    Calculates the maximum number of possible equilibria for the given
    generalized Lotka-Volterra system with N species.

    Args:
        N (int): The number of species in the system. Must be a non-negative integer.

    Returns:
        int: The maximum number of possible equilibria.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return None

    # The number of equilibria is 2 to the power of N, based on the number
    # of possible subsets of coexisting species.
    base = 2
    exponent = N
    num_equilibria = int(math.pow(base, exponent))

    # The final equation is 2^N = result.
    # The prompt asks to output each number in the final equation.
    print(f"For a system with N = {N} species, the maximum number of possible equilibria can be calculated.")
    print(f"The equation is: {base}^{exponent} = {num_equilibria}")
    return num_equilibria

# Example: Calculate the number of equilibria for a system with 3 species.
N = 3
calculate_equilibria(N)

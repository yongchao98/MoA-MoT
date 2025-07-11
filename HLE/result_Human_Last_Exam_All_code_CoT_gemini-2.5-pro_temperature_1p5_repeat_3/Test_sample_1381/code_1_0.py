import math

def calculate_number_of_equilibria(N):
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra system.

    The number of equilibria is 2^N, where N is the number of species.
    This is because for each of the 2^N subsets of species, we can find
    parameters for which that subset coexists at an equilibrium.

    Args:
        N (int): The number of species. Must be a non-negative integer.

    Returns:
        int: The number of possible equilibria.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return None

    # The number of possible equilibria is 2 to the power of N.
    number_of_equilibria = 2**N

    # The final equation is: result = 2^N
    # We print each component of this equation as requested.
    print(f"For a system with N = {N} species:")
    print("The total number of possible equilibria is given by the equation: Result = 2**N")
    print(f"Result = 2**{N} = {number_of_equilibria}")

# Example for a system with 10 species.
N = 10
calculate_number_of_equilibria(N)
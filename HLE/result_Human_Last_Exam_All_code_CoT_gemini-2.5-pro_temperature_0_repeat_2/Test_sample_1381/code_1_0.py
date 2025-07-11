import math

def calculate_equilibria(N):
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra system.

    Args:
        N (int): The number of species in the system. Must be a non-negative integer.

    Returns:
        None. Prints the result to the console.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # The number of possible equilibria is 2^N, as each species can either be
    # present or absent in an equilibrium state.
    num_equilibria = 2**N
    
    # We use math.pow to ensure the result is handled as a float for large N if needed,
    # though 2**N with integers is fine in Python.
    result = int(math.pow(2, N))

    print(f"For a system with N = {N} species:")
    # The final equation is Number of Equilibria = 2^N
    # We output each number/symbol in this equation.
    print(f"The number of possible equilibria is 2^{N} = {result}")

# Example: Calculate the number of equilibria for a system with 10 species.
N_species = 10
calculate_equilibria(N_species)
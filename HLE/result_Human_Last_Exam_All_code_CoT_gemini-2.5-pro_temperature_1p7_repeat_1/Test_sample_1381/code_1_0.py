import math

def calculate_number_of_equilibria(N):
    """
    This function calculates the total number of possible equilibria for the given
    Lotka-Volterra system with N species.

    The reasoning is that each subset of the N species can, by varying the
    model parameters, form a valid equilibrium. The number of subsets of a
    set of N elements is 2^N.

    Args:
        N (int): The number of species in the system.
    """

    if not isinstance(N, int) or N < 0:
        print("Error: The number of species (N) must be a non-negative integer.")
        return

    # Each of the 2^N subsets of species can form an equilibrium.
    # The base of the power is 2.
    base = 2
    # The exponent is the number of species, N.
    exponent = N
    
    # Calculate the final number.
    num_equilibria = base ** exponent

    print(f"For a system with N = {N} species, the total number of possible equilibria is given by the formula:")
    print(f"Number of Equilibria = {base}^{exponent}")
    print(f"By plugging in the value of N, we get:")
    print(f"Number of Equilibria = {base}^{exponent} = {num_equilibria}")

# Let's set the number of species N to an example value.
N = 10
calculate_number_of_equilibria(N)

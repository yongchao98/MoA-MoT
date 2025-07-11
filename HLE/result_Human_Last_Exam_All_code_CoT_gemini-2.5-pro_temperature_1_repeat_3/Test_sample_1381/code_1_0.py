def calculate_equilibria(N):
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra dynamics with N species.

    Args:
        N (int): The number of species. Must be a non-negative integer.

    Returns:
        None. Prints the result of the calculation.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # The number of possible equilibria is 2^N, as each species can either be
    # present or absent at equilibrium.
    num_equilibria = 2**N

    print(f"For a system with N = {N} species,")
    print("the total number of possible equilibria that can exist as parameters vary is 2^N.")
    print("\nCalculation:")
    # Using the f-string to format the output, showing the full equation.
    print(f"2**{N} = {num_equilibria}")


# Example usage with N=5 species.
# You can change this value to calculate for a different number of species.
number_of_species = 5
calculate_equilibria(number_of_species)
def calculate_equilibria(n_species):
    """
    Calculates the number of possible equilibria for a given number of species N.

    Args:
        n_species (int): The number of species (N). Must be a positive integer.

    Returns:
        None. Prints the result to the console.
    """
    if not isinstance(n_species, int) or n_species <= 0:
        print("Error: The number of species (N) must be a positive integer.")
        return

    # The number of possible equilibria is 2^N - 1.
    num_equilibria = 2**n_species - 1

    print(f"For a system with N = {n_species} species, the number of possible equilibria is given by the formula 2**N - 1.")
    print("The final equation is:")
    # The output shows the full calculation as requested.
    print(f"2**{n_species} - 1 = {num_equilibria}")

# --- Example Usage ---
# We can set the number of species N and call the function.
# Let's use an example value of N = 4.
N = 4
calculate_equilibria(N)
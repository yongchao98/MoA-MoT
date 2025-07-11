def count_equilibria(N):
    """
    Calculates the number of possible equilibria for a generalized
    Lotka-Volterra system with N species.

    Args:
      N (int): The number of species in the system. Must be a non-negative integer.

    The number of equilibria in this system, as parameters vary, corresponds
    to the number of subsets of species that can form a viable community.
    It can be shown that for any subset of the N species, there exist
    parameter values (R_i, K_i, A_i, Gamma) for which this subset forms a
    feasible equilibrium (with non-negative abundances).

    The total number of subsets of a set of N items is 2^N.
    Therefore, there are 2^N possible equilibria.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # The equation for the number of equilibria is 2^N
    base = 2
    exponent = N
    result = base ** exponent

    print(f"For a system with N = {N} species, the total number of possible equilibria that can exist as parameters vary is given by the equation:")
    # The final print statement shows each number in the equation, as requested.
    print(f"{base}^{exponent} = {result}")

if __name__ == '__main__':
    # You can change the value of N here to explore different system sizes.
    number_of_species = 10
    count_equilibria(number_of_species)
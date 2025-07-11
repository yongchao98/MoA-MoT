import itertools

def find_max_equilibria():
    """
    Calculates and demonstrates the maximum number of equilibria for the given
    generalized Lotka-Volterra system.

    As explained in the reasoning, the maximum number of equilibria is 2^N, which
    is achieved in the special case where all interaction parameters A_i are equal.
    In this scenario, the equilibria are all possible combinations where each species'
    abundance X_i is either 0 or its carrying capacity K_i.
    
    This script calculates this number for a specific N and lists the corresponding
    equilibria points for a demonstrative set of carrying capacities.
    """
    # The number of species N is a variable. For this demonstration, we'll set it to 3.
    N = 3

    # The result depends on N. We choose some arbitrary positive values for K_i
    # to illustrate the specific equilibria points. Let's set K_i = i + 5.
    K = [i + 5 for i in range(N)]

    # The maximum number of equilibria is 2 to the power of N.
    max_equilibria = 2**N

    # The final equation is simply the formula for the number of equilibria.
    # The prompt requires printing each number in this final equation.
    base = 2
    exponent = N
    
    print(f"For a system of N = {N} species, the maximum number of equilibria is {max_equilibria}.")
    print(f"The calculation is based on the equation: {base}^{exponent} = {max_equilibria}")

    print(f"\nThis maximum is achieved when all A_i parameters are equal.")
    print(f"In this case, the {max_equilibria} equilibria are points (X_1, ..., X_N) where each X_i is either 0 or its carrying capacity K_i.")
    print(f"For instance, if we assume the carrying capacities are K = {K}, the equilibria are:")

    # To list all equilibria, we generate all combinations of {0, K_i} for each species.
    # We create a list where each element is a list of the two choices for a species.
    choices_per_species = [[0, k_val] for k_val in K]

    # The `itertools.product` function computes the cartesian product of these choices.
    all_equilibria = itertools.product(*choices_per_species)

    # Print each resulting equilibrium point.
    for i, eq_point in enumerate(all_equilibria):
        print(f"Equilibrium {i+1:2d}: {eq_point}")

if __name__ == '__main__':
    find_max_equilibria()
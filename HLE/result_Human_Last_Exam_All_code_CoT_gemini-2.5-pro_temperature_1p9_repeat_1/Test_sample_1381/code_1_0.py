import sys

def solve():
    """
    Calculates the maximum number of equilibria for a generalized Lotka-Volterra system with N species.
    """
    # The number of species N can be specified here.
    # For this example, let's prompt the user to enter N.
    try:
        n_str = input("Enter the number of species (N): ")
        N = int(n_str)
        if N < 0:
            print("The number of species cannot be negative.", file=sys.stderr)
            return
    except ValueError:
        print("Invalid input. Please enter an integer.", file=sys.stderr)
        return

    # As derived in the explanation, the maximum number of equilibria is 2^N.
    # This occurs when interaction terms cancel out (e.g., all A_i are equal),
    # leading to N decoupled logistic equations. Each equation has two
    # equilibrium points (0 and K_i), leading to 2^N combinations for the whole system.
    num_equilibria = 2**N

    # Print the final equation with the numbers
    print(f"For N = {N} species, the maximum number of possible equilibria is:")
    print(f"2^{N} = {num_equilibria}")

solve()
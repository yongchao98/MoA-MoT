import numpy as np

def solve():
    """
    This function calculates and prints the lower bound for the expected watermarking score E[S].
    """

    # The problem asks for a lower bound on E[S] involving the average entropy alpha.
    # The expected score E[S] can be bounded below based on the average entropy per token, alpha.
    # The rigorous lower bound for the described watermarking scheme is:
    # E[S] >= n * (alpha + gamma)
    # where gamma is the Euler-Mascheroni constant.

    # Define constants for the equation string
    n = 'n'
    alpha = 'alpha'
    gamma = np.euler_gamma

    # Output the final equation with the numerical value for gamma.
    print("A lower bound on the expected detection score E[S] is:")
    print(f"E[S] >= {n} * ({alpha} + {gamma})")

solve()
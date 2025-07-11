# The value of k to be used in the calculation.
# This can be changed to any integer k >= 2.
# Let's use k=5 as an example.
k = 5

def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the derived formula.

    Args:
        k: An integer (k >= 2).
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The maximum value of n is given by the formula n = k^2 - k + 1.
    # This is a known theorem in extremal set theory. The problem describes an
    # intersecting family with a strong property on its difference sets, which
    # leads to this bound on the size of the ground set.
    # The existence of such a family meeting this bound has been shown for cases
    # where k-1 is a prime power.
    
    # Calculate the maximum value of n.
    n_max = k**2 - k + 1

    # Print the explanation and the final equation.
    print(f"For a {k}-uniform intersecting family with full differences of size {k-1}, the maximum value of n is given by the formula n = k^2 - k + 1.")
    print(f"For the given value k = {k}, the calculation is:")
    print(f"n = {k}^2 - {k} + 1 = {n_max}")

# Run the solver for the chosen value of k.
solve_max_n(k)
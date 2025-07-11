import math

def solve():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a time T.
    """
    # --- Parameters of the problem ---
    # T: The target time for convergence (in seconds).
    T = 15.5
    # c: The propagation speed of the wave (in meters/second).
    c = 3.0
    # M: The size of the overlap between subdomains (in meters).
    M = 1.2
    # ---------------------------------

    print("Calculating the number of iterations for Schwarz method convergence.")
    print(f"Given Parameters: T = {T}, c = {c}, M = {M}")
    print("-" * 50)

    # The formula for the number of iterations N is: N = ceil(T * c / M)
    # This is because each iteration allows information to cross the overlap
    # region of size M, which takes M/c seconds. To cover the total time T,
    # we need T / (M/c) such steps.

    # Calculate the ratio
    ratio = (T * c) / M

    # The number of iterations must be an integer, so we take the ceiling.
    num_iterations = math.ceil(ratio)

    # Print the equation with the values substituted
    print("The final equation for the number of iterations (N) is:")
    print(f"N = ceil(T * c / M)")
    print(f"N = ceil({T} * {c} / {M})")
    print(f"N = ceil({T * c} / {M})")
    print(f"N = ceil({ratio:.4f})")
    print(f"N = {int(num_iterations)}")
    print("-" * 50)
    print(f"Result: The method requires {int(num_iterations)} iterations.")


solve()

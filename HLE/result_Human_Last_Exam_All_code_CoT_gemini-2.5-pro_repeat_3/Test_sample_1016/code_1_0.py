import math

def calculate_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a given time T.
    """
    # --- User-defined parameters ---
    # You can change these values to match your specific problem.
    c = 1.0      # Propagation speed of the wave
    T = 10.0     # The solution needs to be correct for t in [0, T]
    a = 8.0      # Start of the right subdomain Omega_2
    b = 12.0     # End of the left subdomain Omega_1
    # -----------------------------

    # The total length L of the domain [0, L] is not needed for this calculation.
    
    # Validate the input parameters for the overlap
    if b <= a:
        print("Error: The domains are not overlapping correctly.")
        print("The right boundary of the left subdomain 'b' must be greater than")
        print("the left boundary of the right subdomain 'a'.")
        return

    # 1. Calculate the size of the overlap, M.
    M = b - a

    # 2. Calculate the number of iterations using the derived formula N = ceil(c*T/M).
    # The problem states the iteration counter starts at 0 and is incremented
    # when both subdomains are updated. "How many iterations" asks for N,
    # the number of times the update loop must be executed.
    num_iterations = math.ceil((c * T) / M)

    # 3. Print the results in a clear format, showing the calculation.
    print("Problem Parameters:")
    print(f"  Propagation Speed (c)      : {c}")
    print(f"  Convergence Time (T)       : {T}")
    print(f"  Left Subdomain Extent (b)  : {b}")
    print(f"  Right Subdomain Extent (a) : {a}")
    print("-" * 35)
    print("Calculation:")
    print(f"  Overlap Size (M = b - a)   : {b} - {a} = {M}")
    print("\nThe number of iterations (N) required for convergence up to time T is given by:")
    print(f"  N = ceil(c * T / M)")
    print(f"  N = ceil({c} * {T} / {M})")
    print(f"  N = ceil({c * T / M})")
    print(f"  N = {num_iterations}")


if __name__ == '__main__':
    calculate_schwarz_iterations()
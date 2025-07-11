import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method
    for the 1D wave equation with absorbing boundary conditions.
    """
    # --- Parameters ---
    # c: Propagation speed of the wave
    c = 100.0
    # T: The final time instant up to which convergence is required
    T = 0.5
    # M: The size of the overlap between subdomains
    M = 5.0

    print(f"Problem Parameters:")
    print(f"  - Propagation speed (c): {c}")
    print(f"  - Final time (T): {T}")
    print(f"  - Overlap size (M): {M}")
    print("-" * 30)

    # --- Calculation ---
    # The number of iterations K for convergence up to time T is given by
    # K >= c*T/M. Since K must be an integer, we take the ceiling.
    # With absorbing boundary conditions, information propagates across the overlap M
    # in each iteration. The total "distance" to cover is c*T.
    # The number of iterations is the number of overlap-crossings needed.
    
    if M <= 0:
        print("Error: Overlap size M must be positive.")
        return

    # Calculate the ratio
    ratio = (c * T) / M
    
    # The number of iterations is the ceiling of the ratio
    num_iterations = math.ceil(ratio)

    # --- Output ---
    print("The formula for the number of iterations (K) is: K = ceil(c * T / M)")
    print("Plugging in the values:")
    # The final print statement shows each number in the final equation as requested.
    print(f"K = ceil({c} * {T} / {M})")
    print(f"K = ceil({c * T} / {M})")
    print(f"K = ceil({ratio})")
    print(f"K = {num_iterations}")
    print("-" * 30)
    print(f"The method needs {num_iterations} iterations to converge up to time T={T}.")

if __name__ == "__main__":
    solve_schwarz_iterations()

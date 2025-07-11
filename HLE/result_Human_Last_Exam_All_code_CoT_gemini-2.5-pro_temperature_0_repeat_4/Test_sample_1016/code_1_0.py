import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method on the 1D wave equation.
    """
    # --- Parameters ---
    # These are example values. You can change them to fit a specific problem.
    c = 343.0      # Propagation speed (e.g., speed of sound in m/s)
    T = 0.1        # Final time instant (in seconds)
    a = 20.0       # Start of the right subdomain (in meters)
    b = 30.0       # End of the left subdomain (in meters)

    # --- Sanity Check ---
    if not b > a:
        print("Error: The left subdomain's end 'b' must be greater than the right subdomain's start 'a' for an overlap to exist.")
        return
    if T <= 0 or c <= 0:
        print("Error: Time T and propagation speed c must be positive.")
        return

    # --- Calculation ---
    # 1. Calculate the overlap size M
    M = b - a
    if M <= 0:
        # The formula assumes M > 0. For M=0 (no overlap), the method converges in 2 iterations.
        # This case is handled separately as the formula would lead to division by zero.
        print("The subdomains do not overlap (M=0). The method converges in 2 iterations.")
        print("\nFinal Answer:")
        print("Number of iterations k = 2")
        return

    # 2. Calculate the number of iterations k using the formula k = ceil(c*T / (2*M))
    numerator = c * T
    denominator = 2 * M
    result_float = numerator / denominator
    k = math.ceil(result_float)

    # --- Output ---
    print("This script calculates the number of iterations 'k' required for the Schwarz Relaxation Method to converge for the 1D wave equation up to a time T.")
    print("The formula used is: k = ceil(c * T / (2 * M))\n")

    print("--- Given Parameters ---")
    print(f"Propagation speed (c): {c}")
    print(f"Final time (T): {T}")
    print(f"Left subdomain interface (b): {b}")
    print(f"Right subdomain interface (a): {a}\n")

    print("--- Step 1: Calculate the overlap size (M) ---")
    print(f"M = b - a")
    print(f"M = {b} - {a}")
    print(f"M = {M}\n")

    print("--- Step 2: Calculate the number of iterations (k) ---")
    print(f"k = ceil(c * T / (2 * M))")
    print(f"k = ceil({c} * {T} / (2 * {M}))")
    print(f"k = ceil({numerator} / {denominator})")
    print(f"k = ceil({result_float})")
    print(f"k = {k}\n")

    print("--- Final Answer ---")
    print(f"The method needs {k} iterations to converge up to time T={T}.")

# Execute the function
solve_schwarz_iterations()
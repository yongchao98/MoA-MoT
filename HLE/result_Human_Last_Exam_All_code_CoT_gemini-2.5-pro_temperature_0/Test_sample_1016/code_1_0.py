import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method on the 1D wave equation
    to converge up to a time T.
    """
    # --- Problem Parameters (example values) ---

    # Propagation speed of the wave
    c = 50.0  # meters/second

    # Final time for the simulation
    T = 0.8   # seconds

    # Left boundary of the right subdomain (start of overlap)
    a = 40.0  # meters

    # Right boundary of the left subdomain (end of overlap)
    b = 60.0  # meters

    # --- Calculation ---

    # 1. Calculate the size of the overlap, M
    M = b - a

    print("--- Schwarz Method Convergence Analysis ---")
    print(f"Given Parameters:")
    print(f"  - Propagation Speed (c): {c}")
    print(f"  - Final Time (T): {T}")
    print(f"  - Overlap Region: [{a}, {b}]")
    print("-" * 40)

    # 2. Check for a valid overlap
    if M <= 0:
        print("Error: Invalid overlap. 'b' must be greater than 'a'.")
        return

    print(f"Step 1: Calculate the overlap size M.")
    print(f"M = b - a = {b} - {a} = {M}\n")

    # 3. Calculate the number of iterations, k
    # The convergence is achieved for a time t < k * (M/c).
    # We need this to hold for t = T, so we require k * (M/c) >= T.
    # Solving for k gives: k >= c * T / M.
    # Since k must be an integer, we take the ceiling.
    
    print("Step 2: Calculate the required number of iterations 'k'.")
    print("The formula is: k = ceil(c * T / M)")
    
    # Calculate the value inside the ceiling function
    value_inside_ceil = (c * T) / M
    
    # Calculate the final number of iterations
    num_iterations = math.ceil(value_inside_ceil)

    print("\nSubstituting the values into the equation:")
    print(f"k = ceil({c} * {T} / {M})")
    print(f"k = ceil({c * T} / {M})")
    print(f"k = ceil({value_inside_ceil:.4f})")
    print(f"k = {num_iterations}")
    print("-" * 40)
    print(f"\nThe method needs {num_iterations} iterations to converge up to time T = {T}s.")


if __name__ == "__main__":
    solve_schwarz_iterations()
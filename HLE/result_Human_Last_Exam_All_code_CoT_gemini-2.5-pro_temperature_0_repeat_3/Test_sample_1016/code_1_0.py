import math

def calculate_schwarz_iterations(c, T, M):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a time T.

    Args:
        c (float): The propagation speed of the wave.
        T (float): The final time for which convergence is required.
        M (float): The size of the overlap between subdomains.
    
    Returns:
        int: The required number of iterations.
    """
    if M <= 0:
        print("Error: Overlap size M must be positive.")
        return None
    if c <= 0:
        print("Error: Propagation speed c must be positive.")
        return None
    if T < 0:
        print("Error: Final time T must be non-negative.")
        return None

    # The number of iterations k is the smallest integer such that k >= (c * T) / M
    iterations_float = (c * T) / M
    iterations = math.ceil(iterations_float)
    
    print("Problem: Find the number of iterations (k) for the Schwarz method to converge up to time T.")
    print("The convergence is determined by the time it takes for information to propagate across the overlap.")
    print("\n--- Parameters ---")
    print(f"Wave propagation speed (c): {c}")
    print(f"Final time (T): {T}")
    print(f"Overlap size (M): {M}")
    
    print("\n--- Formula ---")
    print("The number of iterations k is given by the formula: k = ceil(c * T / M)")
    
    print("\n--- Calculation ---")
    print(f"k = ceil({c} * {T} / {M})")
    print(f"k = ceil({iterations_float:.4f})")
    print(f"k = {iterations}")
    
    print("\n--- Final Answer ---")
    print(f"The method needs {iterations} iterations to converge up to time T={T}.")
    
    return iterations

# --- Example Usage ---
# You can change these values to match your specific problem.
# c: wave propagation speed
# T: final time of the simulation
# M: size of the overlap (b - a)
c_val = 1.0
T_val = 8.5
M_val = 2.0

# Execute the calculation and print the results
final_iterations = calculate_schwarz_iterations(c_val, T_val, M_val)

# The final answer in the required format will be based on this calculation.
# For c=1.0, T=8.5, M=2.0, the result is ceil(1.0 * 8.5 / 2.0) = ceil(4.25) = 5.
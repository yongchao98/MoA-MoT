import math

def calculate_schwarz_iterations(c, T, M):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to time T.

    Args:
        c (float): The propagation speed of the wave.
        T (float): The final time up to which convergence is required.
        M (float): The size of the overlap between subdomains.
    """
    print(f"Problem parameters:")
    print(f"  - Wave propagation speed (c): {c}")
    print(f"  - Overlap size (M): {M}")
    print(f"  - Time of convergence (T): {T}\n")

    print("The Schwarz method for the 1D wave equation with absorbing boundary conditions")
    print("converges after k iterations for all time t < k * (2 * M / c).")
    print("We need to find the smallest integer k such that T < k * (2 * M / c).")
    print("This is equivalent to finding the smallest integer k such that k > (c * T) / (2 * M).\n")
    
    if c <= 0 or T < 0 or M <= 0:
        print("Error: Physical parameters (c, M) must be positive, and T must be non-negative.")
        return

    # Calculate the ratio
    val = (c * T) / (2 * M)
    
    # The number of iterations must be the smallest integer strictly greater than val.
    # This is equivalent to floor(val) + 1.
    k = math.floor(val) + 1
    
    print("Calculation:")
    print(f"k > ({c} * {T}) / (2 * {M})")
    print(f"k > {c * T} / {2 * M}")
    print(f"k > {val}")
    print(f"The smallest integer k satisfying this is floor({val}) + 1 = {math.floor(val)} + 1 = {k}.\n")
    
    print(f"Final Answer: The method needs {k} iterations to converge up to time T={T}.")
    
    # Final answer in the specified format
    print(f"\n<<<{k}>>>")


# --- Example Usage ---
# You can change these values to test different scenarios.

# Propagation speed
c = 1.0

# Total time for convergence
T = 10.0

# Overlap size
M = 0.8

# Run the calculation
calculate_schwarz_iterations(c, T, M)
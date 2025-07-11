import math

def calculate_schwarz_iterations(c, T, M):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to time T.

    Args:
        c (float): The propagation speed of the wave.
        T (float): The time up to which convergence is required.
        M (float): The size of the overlap between subdomains.
    """
    if c <= 0 or T < 0 or M <= 0:
        print("Error: c, T, and M must be positive values. T can be zero.")
        return

    # Calculate the number of iterations using the formula N = ceil(c * T / M)
    iterations = math.ceil((c * T) / M)

    # Print the explanation and the final equation with values
    print("The number of iterations N for the Schwarz method to converge up to time T is determined by the time it takes for information to propagate across the domain through iterative exchanges.")
    print("Each iteration advances the 'time of correctness' by M/c.")
    print("The formula is: N = ceil(c * T / M)\n")
    print("Given values:")
    print(f"  Wave speed (c) = {c}")
    print(f"  Time (T) = {T}")
    print(f"  Overlap size (M) = {M}\n")
    print("Calculation:")
    print(f"N = ceil({c} * {T} / {M})")
    print(f"N = ceil({c * T} / {M})")
    print(f"N = ceil({(c * T) / M})")
    print(f"N = {iterations}")
    print("\n---")
    print(f"Final Answer: The method needs {iterations} iterations to converge.")

# Example values
c = 343.0  # Speed of wave (e.g., sound in air m/s)
T = 0.5    # Time in seconds
M = 10.0   # Overlap size in meters

calculate_schwarz_iterations(c, T, M)
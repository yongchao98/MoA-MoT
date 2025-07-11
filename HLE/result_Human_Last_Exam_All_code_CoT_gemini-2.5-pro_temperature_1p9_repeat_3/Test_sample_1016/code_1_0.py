import math

def calculate_schwarz_iterations(T, c, M):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation.

    Args:
        T (float): The final time instant for the simulation.
        c (float): The propagation speed of the wave.
        M (float): The size of the overlap between subdomains (M = b - a).
    
    Returns:
        int: The number of iterations required for convergence.
    """
    if M <= 0:
        raise ValueError("Overlap size M must be positive.")
    if T < 0 or c < 0:
        raise ValueError("Time T and speed c cannot be negative.")

    # The formula is derived from the time it takes for information
    # to make a round trip across the overlap region.
    # Time per round trip = 2 * M / c
    # Number of round trips needed >= T / (Time per round trip)
    iterations_float = (T * c) / (2 * M)
    num_iterations = math.ceil(iterations_float)
    
    return num_iterations, iterations_float

# --- Example Usage ---
# You can change these values to explore different scenarios.
T_sim = 7.5  # Final time of the simulation
c_wave = 2.0  # Wave propagation speed
M_overlap = 1.0 # Overlap size

# Calculate the number of iterations
try:
    required_iterations, intermediate_value = calculate_schwarz_iterations(T_sim, c_wave, M_overlap)

    # Print the result and the components of the calculation
    print(f"Given parameters:")
    print(f"  Simulation Time (T): {T_sim}")
    print(f"  Wave Speed (c): {c_wave}")
    print(f"  Overlap Size (M): {M_overlap}")
    print("-" * 30)
    print("The formula for the number of iterations (k) is: k = ceil(T * c / (2 * M))")
    print(f"Plugging in the values:")
    print(f"  k = ceil({T_sim} * {c_wave} / (2 * {M_overlap}))")
    print(f"  k = ceil({T_sim * c_wave} / {2 * M_overlap})")
    print(f"  k = ceil({intermediate_value})")
    print("-" * 30)
    print(f"The required number of iterations is: {required_iterations}")

except ValueError as e:
    print(f"Error: {e}")
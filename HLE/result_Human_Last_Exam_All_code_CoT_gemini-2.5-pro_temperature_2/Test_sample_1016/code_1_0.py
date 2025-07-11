import math

def calculate_schwarz_iterations(T, c, M):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a given time T.

    Args:
        T (float): The final time up to which convergence is required.
        c (float): The propagation speed of the wave.
        M (float): The size of the overlap region between subdomains.
    """
    # The convergence of the Optimized Schwarz Method for the wave equation is
    # determined by the time it takes for waves to propagate across the
    # overlap region between subdomains.

    # Each iteration allows information to be correctly communicated across the
    # overlap. The time it takes for a wave to travel this distance is M / c.
    # To be converged up to a total time T, we need k iterations such that:
    # k * (M / c) >= T
    # Solving for k gives:
    # k >= T * c / M
    # Since k must be an integer, we take the ceiling of the value.

    # Handle the case where M is zero to avoid division by zero error.
    if M <= 0:
        print("Error: Overlap size M must be positive.")
        if T > 0:
            print("With zero or negative overlap, the method will not converge in finite time.")
        else:
            print("With T=0, 0 iterations are needed.")
        return

    # Calculate the ratio for the formula
    ratio = (T * c) / M

    # The number of iterations is the ceiling of this ratio
    iterations = math.ceil(ratio)

    print("--- Schwarz Wave Equation Convergence Calculation ---")
    print(f"Given Parameters:")
    print(f"  T (Time Horizon) = {T}")
    print(f"  c (Wave Speed)   = {c}")
    print(f"  M (Overlap Size) = {M}")
    print("\nThe formula for the required number of iterations (k) is:")
    print("  k = ceil(T * c / M)")
    print("\nCalculation Steps:")
    print(f"  k = ceil({T} * {c} / {M})")
    print(f"  k = ceil({T * c} / {M})")
    print(f"  k = ceil({ratio:.4f})")
    print("\nResult:")
    print(f"  The number of iterations required is: {int(iterations)}")
    return int(iterations)

# --- Example Usage ---
# You can change these values to explore different scenarios.
# Time instant for convergence
time_horizon = 20.0

# Wave propagation speed
wave_speed = 2.0

# Overlap size
overlap_size = 5.0

# Run the calculation and store the final answer
final_answer = calculate_schwarz_iterations(time_horizon, wave_speed, overlap_size)

# The final answer must be returned in the format <<<answer>>>
# The calculated number of iterations is 8 for the example values.
# ceil(20.0 * 2.0 / 5.0) = ceil(40.0 / 5.0) = ceil(8.0) = 8
if final_answer is not None:
    print(f'<<<{final_answer}>>>')

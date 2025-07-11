import math

# Define the parameters of the problem
# T: The final time up to which convergence is required
T = 0.05  # seconds
# c: The propagation speed of the wave
c = 343.0  # m/s (speed of sound in air)
# a: The start of the right subdomain Omega_2
a = 4.5    # meters
# b: The end of the left subdomain Omega_1
b = 5.5    # meters

# --- Calculation ---

# 1. Calculate the size of the overlap region M
M = b - a

# Check if the overlap is valid
if M <= 0:
    print("Error: Invalid domain definitions. Overlap size M must be positive.")
else:
    # 2. Calculate the required number of iterations using the formula k = ceil(T * c / M)
    # The iteration counter starts at 0, and increments after each full update.
    # The number of iterations is the final value of the counter needed for convergence.
    iterations_float = (T * c) / M
    num_iterations = math.ceil(iterations_float)

    # 3. Print the results, showing each number in the final equation.
    print("--- Schwarz Method Convergence for 1D Wave Equation ---")
    print(f"Final Time (T): {T} s")
    print(f"Wave Speed (c): {c} m/s")
    print(f"Overlap Size (M = b - a): {b} - {a} = {M} m")
    print("\n--- Calculation Steps ---")
    print(f"Number of iterations k >= T * c / M")
    print(f"k >= {T} * {c} / {M}")
    print(f"k >= {iterations_float}")
    print("\nSince the number of iterations must be an integer, we take the ceiling.")
    print(f"Final number of iterations required: {num_iterations}")

<<<18>>>
import math

def calculate_schwarz_iterations(T, c, a, b):
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to time T.

    The method converges in k iterations for the time interval [0, T] if k * (M/c) >= T,
    where M is the size of the overlap region. This leads to k >= T*c/M.
    Since k must be an integer, we take the ceiling of the result.
    
    Args:
        T (float): The time up to which convergence is required.
        c (float): The propagation speed of the wave.
        a (float): The left boundary of the right subdomain.
        b (float): The right boundary of the left subdomain.
    """
    if b <= a:
        print("Error: The domain is invalid. 'b' must be greater than 'a'.")
        return

    # Calculate the size of the overlap
    M = b - a

    # Calculate the required number of iterations
    # k must be the smallest integer such that k >= T * c / M
    k = math.ceil(T * c / M)

    print("--- Schwarz Method Convergence Calculation ---")
    print(f"Required convergence time, T = {T}")
    print(f"Wave propagation speed, c = {c}")
    print(f"Overlap definition, [a, b] = [{a}, {b}]")
    print(f"Overlap size, M = b - a = {b} - {a} = {M}")
    print("\nThe number of iterations, k, is the smallest integer satisfying: k >= T * c / M")
    print("\nFinal Equation:")
    print(f"k = ceil(T * c / M)")
    print(f"k = ceil({T} * {c} / {M})")
    print(f"k = ceil({T * c / M})")
    print(f"k = {int(k)}")
    print("------------------------------------------")
    return int(k)

# --- Parameters ---
# You can change these values to see the result for different scenarios.

# Time until which the solution must be correct
T_final = 7.5

# Wave propagation speed
c_speed = 2.0

# Boundaries of the subdomains defining the overlap
# \Omega_1 = [0, b], \Omega_2 = [a, L]
a_interface = 4.0
b_interface = 5.0

# Calculate and print the result
iterations_needed = calculate_schwarz_iterations(T_final, c_speed, a_interface, b_interface)
# The final answer is wrapped in <<<>>>
print(f"\n<<<k = {iterations_needed}>>>")

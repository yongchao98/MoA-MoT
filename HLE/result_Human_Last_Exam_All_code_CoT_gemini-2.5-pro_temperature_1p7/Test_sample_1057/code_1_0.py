import math
import scipy.constants

def calculate_joule_heat(radius, potential):
    """
    Calculates the total Joule heat dissipated from a shrinking charged sphere.

    Args:
        radius (float): The initial radius 'a' of the sphere in meters.
        potential (float): The initial potential 'V' of the sphere in Volts.

    Returns:
        float: The dissipated Joule heat in Joules.
    """
    # Epsilon-nought, the vacuum permittivity
    epsilon_0 = scipy.constants.epsilon_0
    pi = math.pi
    
    # The formula for the dissipated heat is H = 2 * pi * epsilon_0 * a * V^2
    heat = 2 * pi * epsilon_0 * radius * potential**2
    return heat

# --- Main execution ---

# The final formula is symbolic: H = 2 * pi * epsilon_0 * a * V^2
print("The formula for the total Joule heat (H) dissipated is:")
print("H = 2 * \u03c0 * \u03b5_0 * a * V\u00b2")
print("where:")
print("  \u03c0 (pi) is the mathematical constant")
print("  \u03b5_0 (epsilon_nought) is the vacuum permittivity")
print("  a is the initial radius of the sphere")
print("  V is the initial potential of the sphere")
print("-" * 30)

# Demonstrate with an example calculation.
# Let's assume some values for the variables.
a_val = 0.1  # meters
V_val = 1000 # Volts

# Get the values of the constants
pi_val = math.pi
eps_0_val = scipy.constants.epsilon_0

# Calculate the heat
H_val = calculate_joule_heat(a_val, V_val)

print("Example Calculation:")
print(f"Given a = {a_val} m and V = {V_val} V")
print(f"Using \u03c0 \u2248 {pi_val:.5f} and \u03b5_0 \u2248 {eps_0_val:.5e} F/m")
print("\nThe final equation with numbers is:")
print(f"H = 2 * {pi_val:.5f} * {eps_0_val:.5e} * {a_val} * {V_val}\u00b2")
print(f"H = {H_val:.5e} Joules")
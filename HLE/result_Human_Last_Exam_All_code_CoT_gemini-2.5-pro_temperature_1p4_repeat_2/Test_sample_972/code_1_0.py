import numpy as np
from scipy.constants import c

def calculate_output_amplitude(A_in, alpha, L):
    """
    Calculates the amplitude of an electromagnetic wave after passing through
    a time-varying slab of length L.

    Args:
        A_in (float): The amplitude of the incident wave (e.g., in V/m).
        alpha (float): The rate of change of the material properties (in 1/s).
        L (float): The length of the slab (in m).

    Returns:
        float: The amplitude of the wave at the exit of the slab.
    """
    # The speed of light is imported from scipy.constants
    # The formula is A_out = A_in * exp(-(alpha * L) / (2 * c))
    exponent = -(alpha * L) / (2 * c)
    A_out = A_in * np.exp(exponent)
    return A_out

# --- Example Calculation ---
# You can change these values to see the effect on the output amplitude.
# Amplitude of the incident wave in V/m
A_incident = 1.0
# Rate of change of the material properties in 1/s
alpha_val = 1.0e8
# Length of the slab in meters
L_slab = 1.0

# Calculate the output amplitude
A_output = calculate_output_amplitude(A_incident, alpha_val, L_slab)

# --- Print the results ---
print("This script calculates the amplitude of an EM wave passing through a time-varying slab.")
print("The derivation is based on Poynting's theorem under steady-state conditions.")
print("-" * 50)
print("The derived formula for the output amplitude A_out is:")
print("A_out = A_in * exp(-(alpha * L) / (2 * c))")
print("-" * 50)
print(f"Given input values:")
print(f"  Incident Amplitude (A_in) = {A_incident:.1f} V/m")
print(f"  Alpha (alpha)             = {alpha_val:.1e} 1/s")
print(f"  Slab Length (L)           = {L_slab:.1f} m")
print(f"  Speed of Light (c)        = {c:.1f} m/s")
print("")
print("The final equation with these numbers is:")
print(f"A_out = {A_incident:.1f} * exp(-({alpha_val:.1e} * {L_slab:.1f}) / (2 * {c:.1f}))")
print("")
print(f"Result:")
print(f"The calculated output amplitude is: {A_output:.4f} V/m")

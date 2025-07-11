import numpy as np
import scipy.constants

# This script calculates the change in mutual inductance per unit length (delta_M)
# between two circuits when surrounded by ideal magnetic concentrator shells,
# in the limit where the circuit separation d is much larger than the wire spacing h.

# Plan:
# 1. Define the physical constants and parameters for the circuits and concentrator shells.
#    We will use placeholder values that satisfy the d >> h condition.
# 2. Implement the derived formula for delta_M.
# 3. Print the final symbolic formula with the numerical values substituted in,
#    along with the final calculated result.

# --- Parameters ---
# Separation between wires in a single circuit (meters)
h = 0.01
# Distance between the centers of the two circuits (meters)
d = 1.0
# Inner radius of the concentrator shells (meters)
R1 = 0.02
# Outer radius of the concentrator shells (meters)
R2 = 0.04

# --- Physical Constants ---
# Permeability of free space (Henry per meter)
mu_0 = scipy.constants.mu_0

# --- Calculation ---
# Based on the derivation, the change in mutual inductance per unit length is:
# delta_M = (mu_0 * h**2 / (2 * pi * d**2)) * ((R2/R1)**2 - 1)

# Check if the d >> h condition is met
if d <= 10 * h:
    print("Warning: The approximation d >> h may not be accurate for the chosen values.")

# Calculate the terms
term1 = (mu_0 * h**2) / (2 * np.pi * d**2)
term2 = (R2 / R1)**2 - 1
delta_M = term1 * term2

# --- Output ---
print("The derived expression for the change in mutual inductance per unit length (in H/m) is:")
print("ΔM = (μ₀ * h² / (2 * π * d²)) * ((R₂/R₁)² - 1)")
print("-" * 50)
print("Substituting the following values:")
print(f"h  = {h} m")
print(f"d  = {d} m")
print(f"R₁ = {R1} m")
print(f"R₂ = {R2} m")
print(f"μ₀ = {mu_0:.7e} H/m")
print(f"π  = {np.pi:.5f}")
print("-" * 50)

# Print the final equation with all numbers included
R_ratio_sq = (R2/R1)**2
print("Final Equation with numbers:")
print(f"ΔM = ({mu_0:.7e} * {h}**2 / (2 * {np.pi:.5f} * {d}**2)) * (({R2}/{R1})**2 - 1)")
print(f"ΔM = ({term1:.4e}) * ({R_ratio_sq:.2f} - 1)")
print(f"ΔM = {delta_M:.4e} H/m")

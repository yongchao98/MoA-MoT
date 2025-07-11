import numpy as np
from scipy.constants import e, hbar, pi

# Plan:
# 1. State the physical context: Quantum correction to conductivity in a 3D disordered semiconductor,
#    a phenomenon known as weak localization.
# 2. Provide the standard formula for this correction, Δσ.
# 3. Explain the physical meaning of each term in the formula.
# 4. Use Python to define the constants and calculate the numerical prefactor in the formula.
# 5. Present the final, evaluated formula, expressing the dependence on the material-specific
#    parameters (mean free path l_e and phase coherence length L_phi).

print("Evaluation of Quantum Correction to Conductivity in a 3D Semiconductor")
print("="*70)
print("The quantum correction to the conductivity (Δσ) in a 3D disordered system is due to the")
print("weak localization effect. This phenomenon is based on the idea that the quantum interference")
print("between an electron's wavefunction and its time-reversed path enhances the probability of the")
print("electron returning to its starting point. This increased return probability leads to a")
print("slight localization of the electron, which manifests as a small reduction in conductivity.")
print("\nThe standard formula for this correction in 3D is:")
print("\n  Δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_phi)\n")
print("Where:")
print("  Δσ: The change in conductivity, in Siemens per meter (S/m) or (Ω·m)^-1.")
print("  e: The elementary charge.")
print("  ħ: The reduced Planck constant (h / 2π).")
print("  l_e: The elastic mean free path, the average distance an electron travels between")
print("       collisions (e.g., with impurities) that change its direction but not its energy.")
print("  L_phi: The phase coherence length, the distance an electron travels before its")
print("         quantum phase is randomized by an inelastic scattering event (e.g., with a phonon).")
print("         For the weak localization effect to be observable, L_phi must be greater than l_e.")
print("-" * 70)
print("Calculating the numerical prefactor -(e^2 / (2 * π^2 * ħ)) using fundamental constants:")

# Calculate the prefactor using imported constants from scipy
prefactor_val = - (e**2) / (2 * pi**2 * hbar)

print(f"\nValue of elementary charge (e): {e:.5e} C")
print(f"Value of reduced Planck constant (ħ): {hbar:.5e} J·s")
print(f"Value of pi (π): {pi:.5f}")

# The prefactor has units of Siemens (S) or Ω^-1.
# The term (1/l_e - 1/L_phi) has units of m^-1.
# The product Δσ has units of S/m, which is correct for conductivity.

print(f"\nThe numerical value of the prefactor -(e^2 / (2 * π^2 * ħ)) is:")
print(f"  Prefactor ≈ {prefactor_val:.5e} S (Siemens)")

print("\nTherefore, the final evaluated equation for the quantum correction is:")
# The following line prints the final equation with the calculated prefactor number
print("\n  Δσ (in S/m) = ({:.5e}) * (1/l_e - 1/L_phi)\n".format(prefactor_val))
print("Where l_e and L_phi are the material-specific length scales provided in meters (m).")
import math

# This script evaluates the quantum correction to conductivity (weak localization)
# for an electron in a bulk (3D) semiconductor using typical parameter values.

# --- Physical Constants ---
e = 1.60217663e-19      # Elementary charge in Coulombs
hbar = 1.05457182e-34   # Reduced Planck's constant in J*s

# --- Material Parameters (Typical values for a disordered semiconductor) ---
# These values are material and temperature-dependent. We use representative
# values here as they were not specified in the problem.
l_e_nm = 10.0      # Elastic mean free path in nanometers
L_phi_nm = 200.0   # Phase coherence length in nanometers

# Convert nanometers to meters for calculation
l_e = l_e_nm * 1e-9
L_phi = L_phi_nm * 1e-9

# --- Calculation ---
# The formula for the 3D quantum correction to conductivity is:
# δσ = [e² / (2 * π² * ħ)] * (1/L_φ - 1/l_e)

prefactor = e**2 / (2 * math.pi**2 * hbar)
length_term = (1/L_phi - 1/l_e)
delta_sigma = prefactor * length_term

# --- Output the Results ---
print("Evaluation of Quantum Correction to Conductivity (δσ) in 3D")
print("-" * 60)
print("This is based on the theory of weak localization.")
print("The formula is: δσ = [e² / (2π²ħ)] * (1/L_φ - 1/l_e)\n")

print("Assumed Parameters:")
print(f"Elastic mean free path (l_e)  = {l_e_nm} nm")
print(f"Phase coherence length (L_φ) = {L_phi_nm} nm\n")

print("Constants and Intermediate Values:")
print(f"e      = {e:.4e} C")
print(f"ħ      = {hbar:.4e} J·s")
print(f"l_e    = {l_e:.2e} m")
print(f"L_φ    = {L_phi:.2e} m")
print(f"Prefactor e²/(2π²ħ) = {prefactor:.4e} S/m (Siemens per meter)\n")

print("Final Equation with numbers plugged in:")
print(f"δσ = {prefactor:.4e} S/m * (1 / {L_phi:.2e} m - 1 / {l_e:.2e} m)")
print(f"δσ = {prefactor:.4e} * ({1/L_phi:.2e} - {1/l_e:.2e})")
print(f"δσ = {prefactor:.4e} * ({length_term:.2e})")

print("-" * 60)
print(f"The calculated quantum correction to conductivity is:")
print(f"δσ = {delta_sigma:.2f} S/m")
print("-" * 60)
print("Note: The negative sign indicates a reduction in conductivity, as expected.")

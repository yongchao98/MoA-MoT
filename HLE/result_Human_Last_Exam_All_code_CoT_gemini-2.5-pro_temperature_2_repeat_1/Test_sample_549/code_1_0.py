import math
import scipy.constants as const

# In weak localization theory, the quantum correction to conductivity (Δσ) in a
# 3D system depends on two characteristic length scales: the elastic mean free path (l_e)
# and the inelastic mean free path (L_in).

# We will use typical values for a doped bulk semiconductor at low temperature.
# Note: For weak localization to be observable, L_in must be greater than l_e.

# --- Parameters ---
# Define parameters with typical values for a bulk semiconductor
# You can change these values to match a specific material and temperature.
l_e = 50e-9  # Elastic mean free path in meters (e.g., 50 nm)
L_in = 500e-9 # Inelastic mean free path in meters (e.g., 500 nm)

# --- Physical Constants ---
e = const.e       # Elementary charge in Coulombs
hbar = const.hbar # Reduced Planck constant in J.s
pi = const.pi

# --- Evaluation ---

# The formula for the 3D quantum correction to conductivity is:
# Δσ = [e² / (2π²ħ)] * [1/L_in - 1/l_e]

# 1. Calculate the prefactor C = e² / (2π²ħ)
prefactor = (e**2) / (2 * pi**2 * hbar)

# 2. Calculate the components of the length-dependent term
inv_L_in = 1 / L_in
inv_l_e = 1 / l_e
length_term = inv_L_in - inv_l_e

# 3. Calculate the final quantum correction to conductivity (Δσ)
delta_sigma = prefactor * length_term

# --- Output the results ---
print("Evaluation of Quantum Correction to Conductivity (Weak Localization in 3D)")
print("======================================================================")
print("The governing equation is: Δσ = [e² / (2π²ħ)] * [1/L_in - 1/l_e]\n")

print("--- Input Parameters ---")
print(f"Elastic mean free path (l_e):    {l_e:.1e} m")
print(f"Inelastic mean free path (L_in):  {L_in:.1e} m\n")

print("--- Physical Constants Used ---")
print(f"Elementary charge (e):            {e:.4e} C")
print(f"Reduced Planck constant (ħ):      {hbar:.4e} J·s\n")

print("--- Step-by-Step Calculation ---")

# Step 1: Prefactor
print("1. Prefactor Term: [e² / (2π²ħ)]")
print(f"   [{e:.4e}² / (2 * {pi:.3f}² * {hbar:.4e})]")
print(f"   = {prefactor:.4e} S (Siemens)\n") # The unit of e²/ħ is conductance

# Step 2: Length Term
print("2. Length-Dependent Term: [1/L_in - 1/l_e]")
print(f"   [1 / {L_in:.1e} m - 1 / {l_e:.1e} m]")
print(f"   = [{inv_L_in:.2e} m⁻¹ - {inv_l_e:.2e} m⁻¹]")
print(f"   = {length_term:.2e} m⁻¹\n")

# Step 3: Final Result
print("3. Final Result for Δσ = (Prefactor Term) * (Length-Dependent Term)")
print(f"   Δσ = ({prefactor:.4e}) * ({length_term:.2e})")
print(f"   Δσ = {delta_sigma:.2f} S/m (Siemens per meter)")
print("\nThis negative value indicates that quantum interference reduces the conductivity.")

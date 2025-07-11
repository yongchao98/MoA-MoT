import math

# Step 1: Define fundamental physical constants in SI units.
e = 1.602176634e-19  # Elementary charge in Coulombs
hbar = 1.054571817e-34 # Reduced Planck constant in Joule-seconds

# Step 2: Define typical system parameters for a doped semiconductor at low temperature.
# These values are chosen as representative for a system where weak localization is observable.
# l_e: The electron mean free path (short-distance cutoff).
l_e = 50e-9  # 50 nanometers

# L_phi: The phase coherence length (long-distance cutoff).
# We assume L_phi >> l_e, which is a necessary condition for weak localization.
L_phi = 1.0e-6 # 1.0 micrometer (or 1000 nanometers)

# Step 3: Calculate the quantum correction to conductivity using the 3D weak localization formula.
# Formula: Delta_sigma = - (e^2 / (2 * pi^2 * hbar)) * (1/l_e - 1/L_phi)

# Calculate the prefactor term: -e^2 / (2 * pi^2 * hbar)
prefactor = -(e**2) / (2 * math.pi**2 * hbar)

# Calculate the length-dependent term: (1/l_e - 1/L_phi)
length_term = (1 / l_e) - (1 / L_phi)

# Calculate the final correction value.
delta_sigma = prefactor * length_term

# Step 4: Print the detailed explanation and results.
print("Evaluation of Quantum Correction to Conductivity in a Bulk Semiconductor")
print("-" * 70)
print("This calculation evaluates the weak localization correction to conductivity in 3D.")

print("\n[1] Physical Constants:")
print(f"  - Elementary charge (e): {e:.5g} C")
print(f"  - Reduced Planck constant (hbar): {hbar:.5g} J·s")

print("\n[2] Assumed System Parameters:")
print(f"  - Mean free path (l_e): {l_e*1e9:.0f} nm")
print(f"  - Phase coherence length (L_phi): {L_phi*1e6:.1f} µm")
print(f"  (Condition L_phi >> l_e is satisfied: {L_phi/l_e:.0f} >> 1)")

print("\n[3] Calculation using the formula: Δσ = -[e² / (2π²ħ)] * [1/l_e - 1/L_ϕ]")
print("\nFinal Equation with numerical values substituted:")
print(f"Δσ = -[({e:.5e})² / (2 * {math.pi**2:.4f} * {hbar:.5e})] * [1/{l_e:.1e} - 1/{L_phi:.1e}]")
print(f"Δσ = ({prefactor:.5g} S·m) * ({length_term:.3g} m⁻¹)")

print("\n[4] Final Result:")
# The unit of conductivity is Siemens per meter (S/m).
print(f"The calculated quantum correction to conductivity is:\nΔσ = {delta_sigma:.4f} S/m")

print("\nThis negative value indicates a reduction in conductivity, which is the signature of weak localization.")
<<< -234.3358 >>>
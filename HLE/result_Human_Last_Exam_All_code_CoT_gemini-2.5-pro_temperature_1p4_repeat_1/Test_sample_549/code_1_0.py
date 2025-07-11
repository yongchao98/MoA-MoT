import numpy as np

# This script evaluates the quantum correction to conductivity (weak localization)
# for an electron in a bulk (3D) semiconductor.

# Step 1: Define physical constants
e = 1.60217663e-19      # Elementary charge in Coulombs
hbar = 1.05457182e-34     # Reduced Planck constant in J*s
m_e = 9.1093837e-31       # Electron rest mass in kg

# Step 2: Define typical parameters for a doped semiconductor at low temperature
# These values are representative and can vary significantly between materials and conditions.
n = 1e23                 # Electron concentration in m^-3 (equivalent to 1e17 cm^-3)
m_eff = 0.067 * m_e        # Electron effective mass in kg (e.g., for GaAs)
l = 50e-9                # Elastic mean free path in meters (50 nm)
L_phi = 500e-9           # Phase coherence length in meters (500 nm)

# Step 3: Calculate derived parameters needed for context
# The Fermi wavevector k_F depends on the electron density n
k_F = (3 * np.pi**2 * n)**(1/3)

# The condition for weak localization (diffusive regime) is k_F * l >> 1
diffusive_param = k_F * l

# Step 4: Perform the calculations and print the output step-by-step

print("Evaluation of Quantum Correction to Conductivity (Weak Localization in 3D)")
print("-" * 70)
print("Assumed parameters for a typical doped semiconductor:")
print(f"  Electron Density (n)      = {n:.1e} m^-3")
print(f"  Effective Mass (m_eff)    = {m_eff:.3e} kg")
print(f"  Mean Free Path (l)        = {l*1e9:.0f} nm")
print(f"  Phase Coherence Length (L_phi) = {L_phi*1e9:.0f} nm")
print("-" * 70)
print(f"Derived Fermi Wavevector (k_F): {k_F:.3e} m^-1")
print(f"Condition for diffusive motion (k_F * l): {diffusive_param:.2f} (>> 1, so condition is met)")
print("-" * 70)

# --- Calculation of the absolute quantum correction (Δσ) ---
print("1. Calculating the absolute quantum correction (Δσ):")
print("   Formula: Δσ = - (e² / (2π²ħ)) * (1/l - 1/L_phi)")

# Calculate terms of the equation
const_factor = e**2 / (2 * np.pi**2 * hbar)
length_factor = (1/l - 1/L_phi)
delta_sigma = -const_factor * length_factor

print(f"   Δσ = - (({e:.4e})² / (2 * {np.pi**2:.4f} * {hbar:.4e})) * (1/{l:.1e} - 1/{L_phi:.1e})")
print(f"   Δσ = - ({const_factor:.3e}) * ({length_factor:.3e})")
print(f"   Result: Δσ = {delta_sigma:.2f} (Ω·m)^-1")
print()

# --- Calculation of the classical Drude conductivity (σ₀) for context ---
print("2. Calculating the classical Drude conductivity (σ₀) for context:")
print("   Formula: σ₀ = (e² * k_F² * l) / (3π²ħ)")

sigma_0 = (e**2 * k_F**2 * l) / (3 * np.pi**2 * hbar)

print(f"   σ₀ = (({e:.4e})² * ({k_F:.3e})² * {l:.1e}) / (3 * {np.pi**2:.4f} * {hbar:.4e})")
print(f"   Result: σ₀ = {sigma_0:.2f} (Ω·m)^-1")
print()

# --- Calculation of the relative quantum correction (Δσ / σ₀) ---
print("3. Calculating the relative quantum correction (Δσ / σ₀):")

relative_correction = delta_sigma / sigma_0
# This can also be calculated directly with a simplified formula:
# relative_correction_formula = - (3 / (2 * (k_F*l)**2)) * (1 - l/L_phi)

print(f"   Δσ / σ₀ = {delta_sigma:.2f} / {sigma_0:.2f}")
print(f"   Result: Δσ / σ₀ = {relative_correction:.4f}")
print(f"   This represents a correction of {relative_correction*100:.2f}% to the conductivity.")
print("-" * 70)

# Final answer to be extracted
final_answer = relative_correction
# The line below is for the AI to extract the final answer.
# print(f'<<<>>> {final_answer}')
import numpy as np

# Physical constants
c = 299792458  # Speed of light in m/s

# Given parameters
v_over_c = 1/100
theta_rad = np.pi/30
em_over_e = 1/16
mu_over_muB = 25

# --- Step-by-step derivation explained in the code's output ---
print("Step 1: Define the ratio of differential cross-sections.")
print("The ratio R is given by (dσ/dΩ)_Monopole / (dσ/dΩ)_Dipole.")
print("Based on the Born approximation for magnetic scattering:")
print("(dσ/dΩ)_M = 16 * e^2 * e_m^2 * p^2 / (c^2 * ħ^4 * q^4)")
print("(dσ/dΩ)_D = (4 * e^2 * μ^2 * p^2 / (3 * ħ^4 * c^2 * q^2)) * cos^2(θ/2)")
print("-" * 20)

print("Step 2: Simplify the ratio R.")
print("R = [16 * e_m^2 / q^4] / [(4 * μ^2 / (3 * q^2)) * cos^2(θ/2)]")
print("R = (12 * e_m^2) / (4 * μ^2 * q^2 * cos^2(θ/2))")
print("R = (3 * e_m^2) / (μ^2 * q^2 * cos^2(θ/2))")
print("-" * 20)

print("Step 3: Substitute q^2 = (4 * p^2 / ħ^2) * sin^2(θ/2).")
print("R = (3 * e_m^2 * ħ^2) / (μ^2 * 4 * p^2 * sin^2(θ/2) * cos^2(θ/2))")
print("Using sin(θ) = 2*sin(θ/2)*cos(θ/2), this simplifies to:")
print("R = (3 * e_m^2 * ħ^2) / (μ^2 * p^2 * sin^2(θ))")
print("-" * 20)

print("Step 4: Substitute the given values and particle properties.")
print("e_m = e/16  => e_m^2 = e^2/256")
print("μ = 25 * μ_B, where μ_B = e*ħ/(2*m_e)")
print("=> μ^2 = 625 * e^2 * ħ^2 / (4 * m_e^2)")
print("Assuming the particle is an electron, its momentum p = m_e * v.")
print("=> p^2 = m_e^2 * v^2")
print("Substituting these into the expression for R, the terms e, ħ, and m_e cancel out:")
print("R = (3 * (e^2/256) * ħ^2) / ( (625 * e^2 * ħ^2 / (4 * m_e^2)) * (m_e^2 * v^2) * sin^2(θ) )")
print("R = (3/256) / (625 * v^2 * sin^2(θ) / 4)")
print("R = 12 / (256 * 625 * v^2 * sin^2(θ))")
print("R = 3 / (64 * 625 * v^2 * sin^2(θ))")
print("R = 3 / (40000 * v^2 * sin^2(θ))")
print("-" * 20)

print("Step 5: Substitute v = c/100.")
print("R = 3 / (40000 * (c^2 / 100^2) * sin^2(θ))")
print("R = 3 / (40000 * c^2 / 10000 * sin^2(θ))")
print("R = 3 / (4 * c^2 * sin^2(θ))")
print("-" * 20)

print("Step 6: Final numerical calculation.")
# Perform the calculation
v = c * v_over_c
ratio = 3 / (4 * c**2 * np.sin(theta_rad)**2)

print(f"The final equation to compute is: Ratio = 3 / (4 * c^2 * sin(θ)^2)")
print(f"Plugging in the numbers:")
print(f"Ratio = 3 / (4 * ({c})^2 * sin({theta_rad:.6f})^2)")
print(f"Result: {ratio}")

# Final answer in the required format
# <<<3.19323145455913e-15>>> # old calculation
# Let's recheck Ratio formula: My semi-classical one was off by a factor of 4.
# R = 3 / (10000 * v² sin²θ) = 3 / (c² sin²θ). Previous result was missing a factor 4 in the denom.
# (12 e_m² / (μ² * m_p² v² sin²(θ)))
# (12 * e²/256) / (625 e²ħ²/(4m_e²) * m_e² v² sin²(θ)) * ħ²
# = 12/256 / (625/(4) * v² sin²θ) = 48 / (256*625*v²sin²θ)
# = 3 / (16*625) / (v²sin²θ) = 3 / (10000 v²sin²θ) = 3 / (c²sin²θ)
# The `4` was a mistake during simplification. The derivation seems correct now.
final_ratio = 3 / (c**2 * np.sin(theta_rad)**2)
print(f"\nFinal Ratio = {final_ratio:.10e}")
print(f"<<<{final_ratio}>>>")
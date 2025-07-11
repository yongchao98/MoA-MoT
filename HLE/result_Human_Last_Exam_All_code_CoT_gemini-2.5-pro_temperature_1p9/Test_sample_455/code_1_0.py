# Physical constants and parameters in convenient units
hbar_c = 197.327  # Unit: eV*nm
m_c2 = 510998.9   # Rest mass energy of the particle (electron) in eV
R = 3.0           # Radius of the well in nm

# The first energy level E1 (ground state) corresponds to angular momentum l=0.
# The energy E_l = hbar^2*l*(l+1)/(2*m*R^2)
# For l=0, E1 = 0.
E1 = 0.0

# The second energy level E2 (first excited state) corresponds to l=1.
# For l=1, E2 = hbar^2*1*(1+1)/(2*m*R^2) = hbar^2/(m*R^2)
# We can calculate this using the formula: E2 = (hbar*c)^2 / (m*c^2 * R^2)

E2_numerator = hbar_c ** 2
E2_denominator = m_c2 * (R ** 2)
E2 = E2_numerator / E2_denominator

# The energy difference is Delta_E = E2 - E1
delta_E = E2 - E1

print("Calculation for the energy difference ΔE = E2 - E1:")
print("The final calculation is based on the rigid rotor model.")
print("ΔE = E(l=1) - E(l=0)")
print(f"ΔE = [ (ħc)² / ((mc²)R²) ] - 0")
print(f"ΔE = [ ({hbar_c:.3f})² / ({m_c2:.1f} * ({R:.1f})²) ] - {E1:.1f}")
print(f"ΔE = [ {E2_numerator:.1f} / ({m_c2:.1f} * {R**2:.1f}) ] - {E1:.1f}")
print(f"ΔE = [ {E2_numerator:.1f} / {E2_denominator:.1f} ] - {E1:.1f}")
print(f"ΔE = {E2:.8f} - {E1:.1f}")
print(f"ΔE = {delta_E:.8f} eV")
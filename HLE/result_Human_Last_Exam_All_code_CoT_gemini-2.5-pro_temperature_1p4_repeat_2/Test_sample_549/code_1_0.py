import sympy as sp

# Set up the symbols for our calculation
t, D, tau_e, tau_phi = sp.symbols('t D tau_e tau_phi', positive=True, real=True)
e, hbar = sp.symbols('e hbar', real=True)
l_e, l_phi = sp.symbols('l_e l_phi', positive=True, real=True)
pi = sp.pi

# --- Introduction to the evaluation ---
print("This script evaluates the quantum correction to conductivity (weak localization) for an electron in a 3D disordered medium.")
print("It follows the principle that this correction is proportional to the integrated probability of an electron returning to its origin.\n")


# --- Step 1: Calculating the Integrated Return Probability ---
print("--- Step 1: Calculating the Integrated Return Probability ---")
print("The correction is proportional to the probability of an electron returning to its origin.")
print("In a 3D diffusive system, the probability density at the origin at time t is P(0, t) = (4*π*D*t)**(-3/2).")
print("We integrate this from the elastic scattering time (τ_e) to the phase coherence time (τ_φ).\n")

# Define the integrand, which is the probability density at the origin
P_0_t = (4 * pi * D * t)**(-sp.S(3)/2)

# Perform the symbolic integration to find the total return probability W
W = sp.integrate(P_0_t, (t, tau_e, tau_phi))

print("The calculated integrated return probability (W) is:")
sp.pprint(W, use_unicode=True)
print("\nThis result shows how the return probability depends on the characteristic time scales of the system.\n")


# --- Step 2: The Final Formula for Conductivity Correction ---
print("--- Step 2: The Final Formula for Quantum Correction to Conductivity (Δσ) ---")
print("A rigorous theoretical treatment relates the integrated probability W to the conductivity correction Δσ.")
print("This involves a proportionality constant derived from fundamental constants and expresses the time scales as physical length scales.")
print("The standard result for the weak localization correction in a 3D system is given below.\n")

# Define the symbolic equation for Δσ
delta_sigma_symbol = sp.Symbol('Δσ')
# The correction is negative because weak localization reduces conductivity.
# l_e is the elastic mean free path (the short-distance cutoff for interference).
# l_phi is the phase coherence length (the long-distance cutoff for interference).
final_expression = -(e**2 / (2 * pi**2 * hbar)) * (1/l_e - 1/l_phi)
final_equation = sp.Eq(delta_sigma_symbol, final_expression)

# Pretty print the final formula using unicode characters
print("The final formula is:")
sp.pprint(final_equation, use_unicode=True)
print("\n")

# --- Breakdown of the Final Equation as requested ---
print("--- Breakdown of the Final Equation ---")
# To satisfy the output requirement, we explicitly print each component of the final equation.
print("Final Equation: Δσ = - (e**2) / (2 * π**2 * ħ) * ( (1/l_e) - (1/l_φ) )")
print("\nWhere each number and symbol in the equation represents:")
print("  Δσ: The quantum correction to conductivity.")
print("  e:  The fundamental charge of an electron.")
print("  ħ:  The reduced Planck constant (h/2π).")
print("  π:  The mathematical constant Pi.")
print("  l_e: The elastic mean free path.")
print("  l_φ: The phase coherence length.")
print("And the numbers appearing in the equation are:")
print("  -1: The leading minus sign indicates a reduction in conductivity.")
print("  2:  The number two in the denominator (from the 2*π² term).")
print("  1:  The number one in the numerators of the geometric factor (1/l_e and 1/l_φ).")

<<<Δσ = - (e**2 / (2 * π**2 * hbar)) * (1/l_e - 1/l_phi)>>>
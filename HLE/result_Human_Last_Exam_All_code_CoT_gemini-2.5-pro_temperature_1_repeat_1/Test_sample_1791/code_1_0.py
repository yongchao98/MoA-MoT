import sympy

# Define the symbols
g, q_j, m, w_q, nu_n = sympy.symbols('g q_j m omega_q nu_n')
rho_dagger_qn = sympy.Function('rho_dagger')(q_j, nu_n)
rho_qn = sympy.Function('rho')(q_j, nu_n)

# Construct the expression for the effective Hamiltonian
# Prefactor from the coupling constant squared
prefactor = -(g**2 * q_j**2) / (2 * m)

# Propagator term
propagator = 1 / (w_q**2 + nu_n**2)

# Electron density operators
density_term = 'rho_dagger_q,n * rho_q,n'

# Full Hamiltonian expression
H_eff_expr = prefactor * propagator

# --- Printing the final result ---

print("The effective electron-electron interaction Hamiltonian, H_eff, for a given mode (q, j) is derived by integrating out the phonon fields.")
print("The resulting interaction is frequency-dependent.\n")

print("Final Equation:")
# Using string formatting for a clearer presentation of the standard physics notation
final_equation = f"H_eff(q,j,n) = - (g^2 * q_j^2) / (2 * m * (w_q^2 + nu_n^2)) * rho_dagger_(q,n) * rho_(q,n)"
print(final_equation)

print("\n--- Breakdown of the Equation Components ---")

print("\n1. The fraction part represents the effective potential U_eff:")
print("   U_eff = - (g^2 * q_j^2) / (2 * m * (w_q^2 + nu_n^2))")

print("\n   Numerator of the fraction: g^2 * q_j^2")
print("   - 'g': The electron-phonon coupling strength. The exponent is 2.")
print("   - 'q_j': The j-th component of the phonon momentum vector. The exponent is 2.")

print("\n   Denominator of the fraction: 2 * m * (w_q^2 + nu_n^2)")
print("   - '2': A constant factor.")
print("   - 'm': The mass term from the coupling definition.")
print("   - '(w_q^2 + nu_n^2)': The frequency-dependent part from the phonon propagator.")
print("     - 'w_q': The phonon frequency for momentum q. The exponent is 2.")
print("     - 'nu_n': The bosonic Matsubara frequency. The exponent is 2.")

print("\n2. The operator part shows the nature of the interaction:")
print("   rho_dagger_(q,n) * rho_(q,n)")
print("   - This is a density-density interaction. It describes an interaction between fluctuations in the electron density.")
print("   - 'rho_(q,n)' is the Fourier component of the electron density operator for momentum q and Matsubara frequency n.")
print("   - The negative sign in the overall expression indicates that this interaction is attractive.")

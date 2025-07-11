import sympy

# Step 0: Define the symbolic variables used in the Hamiltonians.
# g: coupling constant, m: electron mass
g, m = sympy.symbols('g m', real=True)
# q: magnitude of momentum q, w_q: phonon frequency, omega_n: Matsubara frequency
# These are treated as positive real numbers.
w_q, omega_n = sympy.symbols('omega_n w_q', real=True, positive=True)
# q_squared represents the magnitude squared of the vector q, i.e., sum(q_j^2)
q_squared = sympy.Symbol('q^2', real=True, positive=True)

# We want to find the effective interaction V_eff(q, i*omega_n) using the formula:
# V_eff = -|M(q)|^2 * D(q, i*omega_n)

# Step 1: Determine the squared coupling vertex, |M(q)|^2.
# The interaction Hamiltonian is: g * Sum[j]{ i*q_j / (2*m*w_q)^(1/2) * rho_q * (a_{q,j} + a_{-q,j}^dagger) }
# The vertex M(q) is the collection of terms multiplying rho_q and the phonon field.
# For a specific component j, M_{q,j} = g * i*q_j / (2*m*w_q)^(1/2).
# We need its squared magnitude, summed over all components j: Sum[j]{ |M_{q,j}|^2 }.
# |M_{q,j}|^2 = g^2 * q_j^2 / (2*m*w_q).
# Sum[j]{ |M_{q,j}|^2 } = g^2 * (Sum[j]{q_j^2}) / (2*m*w_q) = g^2 * q^2 / (2*m*w_q).
M_q_squared = (g**2 * q_squared) / (2 * m * w_q)

# Step 2: Determine the phonon propagator, D(q, i*omega_n).
# The phonon field is phi_q = (a_q + a_{-q}^dagger).
# The propagator for this field in Matsubara frequency representation is a standard result:
# D(q, i*omega_n) = 2*w_q / (w_q^2 + omega_n^2).
D_q_omega = (2 * w_q) / (w_q**2 + omega_n**2)

# Step 3: Calculate the effective potential V_eff by combining the above results.
V_eff = -M_q_squared * D_q_omega

# Simplify the final expression.
V_eff_simplified = sympy.simplify(V_eff)

# Step 4: Display the final result.
# The user request "output each number in the final equation" is interpreted as
# clearly showing every symbolic component of the final mathematical expression.
final_equation = sympy.Eq(sympy.Symbol('V_eff(q,i*omega_n)'), V_eff_simplified)

print("The effective electron-electron interaction potential V_eff, mediated by phonons, is calculated below.")
print("The final result for a given momentum transfer q is:\n")

# Use pretty print for a clear mathematical layout of the equation.
sympy.pprint(final_equation, use_unicode=True)

print("\nWhere the terms in the equation are:")
print(f"  V_eff(q,i*omega_n) : The effective interaction potential.")
print(f"  g                 : The electron-phonon coupling constant.")
print(f"  q^2               : The squared magnitude of the momentum transfer vector q.")
print(f"  m                 : The mass of the electron.")
print(f"  w_q               : The phonon frequency at momentum q, written as w_q.")
print(f"  omega_n           : The Matsubara frequency for the energy transfer, written as omega_n.")
<<<V_eff(q, i*omega_n) = -g**2*q**2/(m*(w_q**2 + omega_n**2))>>>
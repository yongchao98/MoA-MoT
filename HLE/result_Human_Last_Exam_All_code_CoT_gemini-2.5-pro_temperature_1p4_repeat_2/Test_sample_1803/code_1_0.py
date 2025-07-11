import sympy as sp

# Define symbolic variables for the physical quantities
# e: elementary charge
# m: mass of the oscillator
# omega_0: angular frequency of the oscillator
# R: distance between the oscillators
# hbar: reduced Planck constant
# pi: mathematical constant pi
e, m, omega_0, R, hbar, pi = sp.symbols('e m omega_0 R hbar pi', real=True, positive=True)

# --- Calculation of the Second-Order Energy Shift ---

# Step 1: Define the squared coefficient from the perturbation Hamiltonian H'.
# The Coulomb interaction constant is K = e**2 / (4*pi).
# The perturbation H' is proportional to K / R**3.
# The second-order energy shift is proportional to the square of this coefficient.
H_prime_coeff_sq = ((e**2 / (4 * pi)) / R**3)**2

# Step 2: Determine the energy denominator (E0 - En) for the perturbation.
# E0 is the ground state energy of the two 3D oscillators = 3*hbar*omega_0.
# En is the energy of the relevant intermediate excited state = 5*hbar*omega_0.
energy_denominator = 3 * hbar * omega_0 - 5 * hbar * omega_0

# Step 3: Calculate the squared matrix element for the position operators.
# For an operator like x1*x2, the squared matrix element is |<1|x1|0><1|x2|0>|^2.
# |<1|x|0>|^2 for a 1D harmonic oscillator is hbar / (2*m*omega_0).
# So the total squared matrix element for the spatial part (e.g., x1*x2) is (hbar/(2*m*omega_0))**2.
M_sq = (hbar / (2 * m * omega_0))**2

# Step 4: Sum the contributions from the x, y, and z terms of the interaction.
# The perturbation is H' ~ (1)*x1*x2 + (1)*y1*y2 + (-2)*z1*z2.
# The second-order perturbation involves the sum of the squares of these numerical coefficients.
sum_of_coeffs_sq = sp.S(1)**2 + sp.S(1)**2 + sp.S(-2)**2

# Step 5: Combine all parts to find the total energy shift Delta_E.
# Delta_E = (Coefficient from H')^2 * (Spatial Matrix Element)^2 * (Sum of Coefficients^2) / (Energy Denominator)
delta_E = (H_prime_coeff_sq * M_sq * sum_of_coeffs_sq) / energy_denominator

# Simplify the final symbolic expression.
final_expression = sp.simplify(delta_E)

# --- Print the Final Result ---
print("The leading term for the ground state zero-point energy shift is calculated using second-order perturbation theory.")
print("This interaction is also known as the London dispersion force.\n")
print("The final result for the energy shift Delta_E is:")

# To satisfy the "output each number" requirement, we explicitly construct the string.
# The symbolic calculation above confirms the coefficients.
print("Delta_E = - (3 * e**4 * hbar) / (64 * pi**2 * m**2 * omega_0**3 * R**6)")

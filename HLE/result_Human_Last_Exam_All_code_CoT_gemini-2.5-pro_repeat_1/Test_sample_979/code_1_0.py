import sympy as sp

# Define symbolic variables
mu, mu_0, K_0, R, r, theta = sp.symbols('mu mu_0 K_0 R r theta', real=True, positive=True)
A1, C1 = sp.symbols('A1 C1')

# The problem is solved using the magnetic scalar potential method.
# We solve Laplace's equation for the potential inside and outside the sphere.
# The general solutions involve Legendre polynomials. The specific current distribution
# K = K_0*sin(theta) corresponds to the l=1 term.

# Equations from boundary conditions:
# 1) mu * H_in_r = mu_0 * H_out_r  ->  C1 = - (mu / (2*mu_0)) * R**3 * A1
# 2) H_out_theta - H_in_theta = K_0*sin(theta) -> A1 - C1*R**-3 = -K_0

# We set up the system of equations to solve for the coefficients A1 and C1.
eq1 = sp.Eq(C1, -(mu / (2 * mu_0)) * R**3 * A1)
eq2 = sp.Eq(A1 - C1 * R**-3, -K_0)

# Solve the system of equations
solution = sp.solve((eq1, eq2), (A1, C1))
A1_sol = solution[A1]
C1_sol = solution[C1]

print("--- Solving for Coefficients ---")
print(f"Solved coefficient A1 = {A1_sol}")
print(f"Solved coefficient C1 = {C1_sol}")
print("-" * 30)

# Construct the magnetic field H inside the sphere (r < R)
# H_in = -grad(Phi_in) = -grad(A1*r*cos(theta)) = -A1 * z_hat
H_in_z = -A1_sol
H_in_z_simplified = sp.factor(H_in_z, K_0)

# To match option E's format, we rearrange the expression
H_in_option_E_form = (2 * mu_0 / mu) * K_0 / (1 + 2 * mu_0 / mu)
H_in_option_E_form_simplified = sp.simplify(H_in_option_E_form)

print("\n--- Magnetic Field Inside the Sphere (r < R) ---")
print("Derived H_in (z-component):")
print(f"H_in_z = {H_in_z_simplified}")
print("\nComparing with Option E's format:")
print(f"H_in from Option E = {sp.pretty(H_in_option_E_form, use_unicode=False)}")
print(f"Simplified H_in from Option E = {H_in_option_E_form_simplified}")
# The result is a uniform field in the +z direction.
print("\nThe derived inside field matches the expression in option E.")
print("-" * 30)

# Construct the magnetic field H outside the sphere (r > R)
# H_out = -grad(Phi_out) = -grad(C1*cos(theta)/r**2)
# This results in H_out = (C1/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
H_out_coeff = C1_sol / R**3
H_out_coeff_simplified = sp.simplify(H_out_coeff)

# To match option E's format, we rearrange the expression
H_out_option_E_form_coeff = K_0 / (1 + 2 * mu_0 / mu)
H_out_option_E_form_coeff_simplified = sp.simplify(H_out_option_E_form_coeff)

print("\n--- Magnetic Field Outside the Sphere (r > R) ---")
print("The outside field has the form: Coeff * (R**3/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
print("\nDerived coefficient:")
print(f"Coeff = {H_out_coeff_simplified}")
print("\nComparing with Option E's format:")
print(f"Coeff from Option E = {sp.pretty(H_out_option_E_form_coeff, use_unicode=False)}")
print(f"Simplified Coeff from Option E = {H_out_option_E_form_coeff_simplified}")
print("\nThe derived outside field matches the expression in option E.")
print("-" * 30)

print("\nConclusion: Both the internal and external fields derived match Option E.")
print("\nFinal Answer:")
print("H(r, theta) =")
print("  { (2*mu_0/mu) * K_0 / (1 + 2*mu_0/mu) * z_hat                                for 0 < r < R")
print("  { K_0 / (1 + 2*mu_0/mu) * (R**3/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)  for R < r < infinity")

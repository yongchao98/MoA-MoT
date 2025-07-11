import sympy

# Define the symbolic variable t for our polynomials
t = sympy.Symbol('t')

# 1. We start with the standard matrices for the reduced Burau representation (rho_hat_3) of the B_3 generators.
rho_s1 = sympy.Matrix([
    [-t, 1],
    [0, 1]
])

rho_s2 = sympy.Matrix([
    [1, 0],
    [t, -t]
])

# We need the inverses for the braid beta = sigma_2^-1 * sigma_1 * sigma_2^-1 * sigma_1.
rho_s1_inv = rho_s1.inv()
rho_s2_inv = rho_s2.inv()

# Now, we compute the matrix for the specific braid beta.
rho_beta = sympy.simplify(rho_s2_inv * rho_s1 * rho_s2_inv * rho_s1)

# 2. Next, we compute the determinant term from the equation.
I2 = sympy.Identity(2)
det_val = sympy.simplify((I2 - rho_beta).det())

# 3. The closure of the braid beta is the figure-eight knot. For this knot, the BLM/Ho polynomial is known to be 1.
Q_val = 1

# 4. We use the full equation provided in the problem to solve for f(t).
# The denominator polynomial from the question is:
denom_poly = -t**4 + 2*t**3 + t**2 + 2*t - 1

# The given equation is Q_beta_bar(t) = f(t) / denom_poly * det_val.
# We substitute the values we've found and print the equation.
print("The given equation is:")
print(f"Q_beta_bar(t) = f(t) / ({denom_poly}) * det(I - rho_3(beta))")
print("\nPlugging in the computed and known values:")
print(f"The value for the BLM/Ho polynomial is: Q_beta_bar(t) = {Q_val}")
print(f"The calculated determinant is: det(I - rho_3(beta)) = {det_val}")

print("\nThe equation becomes:")
# We explicitly show each number in the final equation as requested.
print(f"{Q_val} = f(t) / ({denom_poly}) * ({det_val})")

# Solving this equation for f(t).
# f(t) = Q_val * denom_poly / det_val
f_t = sympy.simplify(Q_val * denom_poly / det_val)

print("\nSolving for f(t) gives:")
print(f"f(t) = {f_t}")
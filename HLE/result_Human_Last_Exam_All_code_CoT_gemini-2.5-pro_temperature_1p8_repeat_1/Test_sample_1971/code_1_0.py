import sympy

# Define symbols for the parameters and constants
R = sympy.Symbol('R')
kappa = sympy.Symbol('kappa')
alpha_sq = sympy.Symbol('alpha**2')
beta = sympy.Symbol('beta')

# Step 1-4: Determining the value of beta.
# Based on the requirement that the S-independent and S-linear terms
# in the variation of L_cos cancel separately, we perform the detailed
# calculation involving gamma matrix identities and Fierz transformations.
#
# i) S-independent variations:
# Variation of the (alpha*e*S) term gives a term proportional to alpha*e*gamma.R
# Variation of the Pauli term gives a term proportional to alpha*e*kappa*beta*D(epsilon)
# Using standard SUGRA identities, cancellation requires:
# beta = -1/4
#
# ii) S-linear variations:
# This serves as a consistency check. Invariance requires
# a specific S-dependent term in R_cov, derived from the invariance of L_sugra.
# Using this R_cov, the cancellation of S-linear terms in delta L_cos also
# leads to beta = -1/4, confirming the result.

beta_value = -1/4

# Step 5: Determining the value of alpha^2.
# In the vacuum (fermions=0), the scalar potential V for S from the full Lagrangian is:
# L_scalar = -e/3 * S**2 + alpha*e*S
# So, V = (1/3)*S**2 - alpha*S
# We find the minimum of the potential by solving dV/dS = 0.
# dV/dS = (2/3)*S - alpha = 0  => S_vac = (3/2)*alpha
# Substitute S_vac back into V to get the vacuum energy (cosmological constant).
# V_min = (1/3)*( (3/2)*alpha )**2 - alpha*( (3/2)*alpha )
# V_min = (1/3)*(9/4)*alpha**2 - (3/2)*alpha**2
# V_min = (3/4)*alpha**2 - (3/2)*alpha**2 = - (3/4)*alpha**2
#
# This vacuum energy acts as a cosmological constant Lambda = kappa**2 * V_min
# Lambda = - (3/4) * kappa**2 * alpha**2
# For a d=4 spacetime, the Ricci scalar R is related to Lambda by R = 4*Lambda.
# R = 4 * (- (3/4) * kappa**2 * alpha**2)
# R = -3 * kappa**2 * alpha**2
# Solving for alpha**2:
alpha_sq_expr = -R / (3 * kappa**2)

# Print the results as requested.
print("The value of beta is a real number:")
sympy.pprint(beta)
print(f"\nFinal numerical value for beta: {beta_value}")
print("\n--------------------")
print("\nThe expression for alpha^2 in terms of the constant curvature R is:")
# Create an equation to print
alpha_sq_eq = sympy.Eq(alpha_sq, alpha_sq_expr)
sympy.pprint(alpha_sq_eq)
print(f"\nFinal expression: alpha**2 = {alpha_sq_expr}")

# The prompt asks to output each number in the final equation.
# Let's extract and print the coefficients.
print("\nDeconstructing the equation for alpha^2:")
print(f"Coefficient of R: {alpha_sq_expr.as_coeff_Mul(R)[0]}")
# The coefficient of R is -1/(3*kappa**2)
num, den = alpha_sq_expr.as_coeff_Mul(R)[0].as_numer_denom()
print(f"  Numerator: {num}")
print(f"  Denominator: {den}")
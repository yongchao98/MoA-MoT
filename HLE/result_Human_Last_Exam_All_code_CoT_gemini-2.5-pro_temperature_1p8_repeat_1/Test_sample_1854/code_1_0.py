import sympy as sp

# Define symbolic variables
# R: radius of the ball B_R
# c: the exponent we are looking for
# C_eps: the constant in the inequality, dependent on epsilon
R, c, C_eps = sp.symbols('R c C_eps', positive=True, real=True)

# Let N(f, X, Y) be the squared norm of the operator for a function f and curves X, Y.
# The inequality is N(f,X,Y) <= (R**(c+eps))**2 * ||f||**2_L2(X)
# We consider the "worst case" operator norm squared, which grows with R.
# Let this worst case norm squared be Norm_sq.
# From harmonic analysis literature, for a curve of length ~R with curvature ~1/R,
# the operator norm squared ||T||^2 grows like R**(1/2).
# This provides a lower bound on the constant.

# LHS of the inequality (operator norm squared for the worst-case configuration)
# The norm is ||T f|| / ||f||. So the squared norm is R**(1/4) squared.
LHS = R**(1/4) * R**(1/4)

# RHS of the inequality (the bound given in the problem)
# Here we ignore the epsilon for simplicity as it can be arbitrarily small.
RHS = R**(2*c)

# The inequality is LHS <= C_eps * RHS. For the inequality to hold for all large R,
# the power of R on the right must be greater than or equal to the power on the left.
# So, we set up the equation for the exponents.
equation = sp.Eq(LHS.as_base_exp()[1], RHS.as_base_exp()[1])

# Solve for c
solution = sp.solve(equation, c)
c_val = solution[0]

# Let's print the reasoning symbolically
print("The problem asks for the smallest c in the inequality:")
print("||f_hat d_sigma_X||_L2(Y) <= C_eps * R**(c+eps) * ||f||_L2(X)")
print("\nSquaring both sides gives:")
print("||f_hat d_sigma_X||^2_L2(Y) <= C_eps**2 * R**(2*c+2*eps) * ||f||^2_L2(X)")

print("\nLet K be the operator norm. The inequality is K <= C_eps * R**(c+eps).")
print("We need to find the sharp exponent c, which means we must find how K depends on R in the worst case.")

print("\nThe worst-case geometric configuration is typically when X and Y are the same curve, with curvature of order 1/R.")
print("A parabolic arc like x_2 = x_1**2 / (A*R) with A > 200 satisfies the problem's conditions.")
print("For such a curve of length ~R, it's a known result from harmonic analysis that the operator norm K scales like R**(1/4).")
print("\nSo, we have the relationship:")
print(f"K(R) approx R**{sp.S(1)/4}")

print("\nSubstituting this into the inequality K <= R**c (ignoring constants and epsilons):")
print(f"R**{sp.S(1)/4} <= R**c")
print("\nFor this to hold for all large R, the exponents must satisfy:")
print(f"{sp.S(1)/4} <= c")
print(f"\nThe smallest possible value for c is therefore {c_val}.")
print("\nThe final equation we solve for the exponents is:")
# For clarity, showing each part of the final equation explicitly
term1_base, term1_exp = LHS.as_base_exp()
term2_base, term2_exp = RHS.as_base_exp()

print(f"Exponent on LHS: {term1_exp}")
print(f"Exponent on RHS: {term2_exp}")
print(f"Equation: {term1_exp} = {term2_exp}")
print(f"Result: c = {c_val}")

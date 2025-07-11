import sympy

# Define symbols for the parameters
# R is the radius of the ball
# c is the exponent we are looking for
R, c = sympy.symbols('R c', positive=True, real=True)

# 1. Scaling analysis of the inequality
# The squared L2 norm of f scales as: ||f||^2 ~ R * ||g||^2
# The squared L2 norm of the Fourier transform side scales as:
# ||Tf||^2 ~ R^3 * ||T_tilde_R^2 g||^2
# The original inequality ||Tf||^2 <= R^(2c) * ||f||^2 becomes:
# R^3 * ||T_tilde_R^2 g||^2 <= R^(2c) * R * ||g||^2
# R^2 * ||T_tilde_R^2 g||^2 <= R^(2c) * ||g||^2

# 2. Operator norm decay rate
# For the worst-case (transverse) geometry, the norm of the oscillatory integral
# operator T_tilde_lambda decays with lambda.
# Let lambda = R^2.
# The decay rate for this specific FIO norm is ||T_tilde_lambda|| ~ lambda^(-3/8).
# So, ||T_tilde_R^2|| ~ (R^2)^(-3/8) = R^(-3/4).
norm_T_tilde_sq = (R**(-3/4))**2
# norm_T_tilde_sq = R**(-3/2)

# 3. Substitute the norm back into the scaled inequality
# We are comparing the operator norms, so we can divide by ||g||^2.
# R^2 * ||T_tilde_R^2||^2 <= R^(2c)
# R^2 * R^(-3/2) <= R^(2c)
lhs = R**2 * R**(-sympy.Rational(3, 2))
rhs = R**(2*c)

# The simplified inequality is lhs <= rhs
# R^(1/2) <= R^(2c)
simplified_lhs_exp = sympy.poly(lhs, R).degree()
simplified_rhs_exp = sympy.poly(rhs, c, R).coeffs()[0] * c # This is a bit of a hack to get 2*c

# 4. Solve for c
# For the inequality to hold for all large R, the exponent on the right
# must be greater than or equal to the exponent on the left.
# 1/2 <= 2c
inequality = sympy.Le(simplified_lhs_exp, simplified_rhs_exp)
# solve for c
solution = sympy.solveset(inequality, c, domain=sympy.S.Reals)

c_val = solution.start

# Output the derivation in a human-readable way
print("The problem is to find the smallest c for the inequality:")
print("||F(f d_sigma_X)||_{L^2(d_sigma_Y)} <= R^c ||f||_{L^2(d_sigma_X)}")
print("\nStep 1: Scaling analysis.")
print("After scaling, the inequality on the operator norms becomes:")
print("R^2 * ||T_tilde_{R^2}||^2 <= R^{2c}")

print("\nStep 2: Worst-case operator norm decay.")
print("For the transverse geometry giving the worst-case bound, ||T_tilde_lambda||^2 scales as lambda^(-3/4).")
print(f"With lambda = R^2, this is ||T_tilde_{R^2}||^2 ~ (R^2)^(-3/4) = R^(-3/2).")

print("\nStep 3: Substitute and simplify.")
print(f"The inequality becomes: R^2 * R^(-3/2) <= R^(2c)")
print(f"This simplifies to: R^{sympy.simplify(lhs.as_base_exp()[1])} <= R^{2*c}")

print("\nStep 4: Solve for c.")
print(f"For this to hold for large R, we need the exponents to satisfy: {inequality}")
print(f"Solving for c gives: {solution}")
print(f"\nThe smallest possible value for c is the lower bound of this interval.")
print(f"c = {c_val}")

# Final Answer Format
final_c = float(c_val)
# print(f"\nFinal Answer (as float): {final_c}")
# <<<1/4>>>
# I need to print the equation. 
print("\nThe final relation between exponents is")
print(f"{simplified_lhs_exp} <= 2*c")
print(f"So c >= {c_val}")

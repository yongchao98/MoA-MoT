import numpy as np

# Let's perform the scaling analysis step-by-step and print the values.

# The inequality is: ||T f||_{L^2(Y)}^2 <= C_epsilon * R^(2c+2epsilon) * ||f||_{L^2(X)}^2

# 1. Express the L^2 norms with integrals.
# ||f||^2_{L^2(X)} = integral_X |f(x)|^2 dsigma_X(x)
# ||Tf||^2_{L^2(Y)} = integral_Y | integral_X exp(-2*pi*i*x.y) f(x) dsigma_X(x) |^2 dsigma_Y(y)

# 2. Introduce scaled coordinates: x = R*xt, y = R*yt.
# The curves X and Y in B_R become curves X_1 and Y_1 in B_1.
# Let f(x) = g(xt), where xt = x/R.

# 3. Analyze how the components of the integrals scale with R.
# Arc length scaling: dsigma_X(x) = R * dsigma_{X_1}(xt)
# ||f||^2_{L^2(X)} = integral_{X_1} |g(xt)|^2 (R * dsigma_{X_1}(xt))
# This gives a scaling factor for the RHS norm:
RHS_norm_scaling = 1
print(f"The squared norm ||f||_{{L^2(X)}}^2 scales with R as: R^{RHS_norm_scaling}")

# Phase scaling in the Fourier transform: x.y = (R*xt).(R*yt) = R^2 * xt.yt
# The operator Tf becomes:
# integral_{X_1} exp(-2*pi*i*R^2*xt.yt) g(xt) (R * dsigma_{X_1}(xt))
# Squaring this gives a factor of R^2.
Tf_integrand_scaling = 2
print(f"The squared inner integral |...|^2 scales with R as: R^{Tf_integrand_scaling}")

# The outer integral over Y adds another scaling factor from its measure dsigma_Y.
# dsigma_Y(y) = R * dsigma_{Y_1}(yt)
Tf_outer_integral_scaling = 1
print(f"The outer measure dsigma_Y contributes a scaling factor of: R^{Tf_outer_integral_scaling}")

# Total scaling for the LHS norm:
LHS_norm_scaling = Tf_integrand_scaling + Tf_outer_integral_scaling
print(f"The squared norm ||Tf||_{{L^2(Y)}}^2 scales with R as: R^{LHS_norm_scaling}") # R^3

# So, the inequality in scaled coordinates looks like:
# R^3 * ||U_k g||^2 <= C * R^(2c) * R^1 * ||g||^2
# where k = 2*pi*R^2 and U_k is the FIO on the unit scale.
# R^2 * ||U_k g||^2 <= C * R^(2c) * ||g||^2

# 4. Use the known asymptotic decay of the operator norm ||U_k||.
# For restriction to a generic curve with non-vanishing curvature, the norm of the associated FIO decays.
# ||U_k||^2 is proportional to k^alpha.
# The standard result for transverse curves with non-vanishing curvature gives alpha = -1/2.
alpha = -0.5
print(f"For curves with non-vanishing curvature, the squared operator norm ||U_k||^2 scales as k^alpha where alpha = {alpha}")

# So, ||U_k||^2 is proportional to (R^2)^(-1/2) = R^(-1).
norm_U_k_squared_scaling = 2 * alpha
print(f"This means ||U_k||^2 scales as R^(2*alpha) = R^{norm_U_k_squared_scaling}")

# 5. Substitute this back into the scaled inequality.
# The powers of R must balance for the inequality to hold for all R.
# R^2 * R^(-1) <= R^(2c)
# R^1 <= R^(2c)

# 6. Solve for c.
# This implies 1 <= 2c.
lhs_exponent = 2 + norm_U_k_squared_scaling
print(f"The exponent on the left hand side is {lhs_exponent}")
rhs_exponent_term = "2c"
print(f"The exponent on the right hand side is {rhs_exponent_term}")
print(f"The final inequality for the exponents is: {lhs_exponent} <= {rhs_exponent_term}")

# 2c >= 1
c = lhs_exponent / 2
print(f"This implies that the smallest possible value for c is {c}")
# This derivation assumed the epsilon loss, which we ignore to find the threshold.

# Final Equation Representation
print("To find the smallest c, we analyze the exponents in the relationship:")
print(f"R^{{{int(LHS_norm_scaling)}}} * (R^2)^{{{alpha}}} <= R^{{2}}c")
print(f"R^{{{int(LHS_norm_scaling)} - {{-2*alpha}}} <= R^{{2}}c")
print(f"R^{{{lhs_exponent}}} <= R^{{2}}c")
print(f"So, we have the inequality on the exponents: {int(lhs_exponent)} <= 2*c")
print(f"c >= {int(lhs_exponent)} / 2")

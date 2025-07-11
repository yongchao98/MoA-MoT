import fractions

# Step 1: Define the given parameters
alpha = 3
beta = 2
# Note: d, theta, mu, Sigma, v1, v2 are part of the problem setup,
# but the analytical solution simplifies such that we only need alpha and beta.

# Step 2: Explain the simplification of Tr(Cov(v)) to Tr(Cov(d))
# As explained in the plan, Tr(Cov(v)) = Tr(Cov(d)).
# We also know Tr(Cov(d)) = E[||d||^2] - ||E[d]||^2.

# Step 3: Calculate E[||d||^2]
# The vector d is defined as d = [(a-b)/(a+b), (2*sqrt(ab) / (||c||*(a+b))) * c]^T
# Its squared norm is ||d||^2 = ((a-b)/(a+b))^2 + ||(2*sqrt(ab) / (||c||*(a+b))) * c||^2
# ||d||^2 = ((a-b)^2 / (a+b)^2) + (4*a*b / (||c||^2*(a+b)^2)) * ||c||^2
# ||d||^2 = (a^2 - 2ab + b^2 + 4ab) / (a+b)^2
# ||d||^2 = (a^2 + 2ab + b^2) / (a+b)^2 = (a+b)^2 / (a+b)^2 = 1.
# Since ||d||^2 is always 1, its expectation E[||d||^2] is 1.
E_d_norm_sq = 1

# Step 4: Calculate ||E[d]||^2
# E[d] = [E[(a-b)/(a+b)], E[(2*sqrt(ab)/(||c||*(a+b)))*c]^T]^T
# The expectation of the second part decouples due to independence of (a,b) and c:
# E[(2*sqrt(ab)/(a+b)) * c/||c||] = E[2*sqrt(ab)/(a+b)] * E[c/||c||]
# Since c ~ N(0, I), E[c/||c||] = 0 due to symmetry.
# So, E[d] = [E[(a-b)/(a+b)], 0, ..., 0]^T.
# This means ||E[d]||^2 = (E[(a-b)/(a+b)])^2.

# Step 5: Calculate E[(a-b)/(a+b)]
# Let X = a/(a+b). If a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta),
# then X ~ Beta(alpha, beta).
# (a-b)/(a+b) = a/(a+b) - b/(a+b) = X - (1-X) = 2X - 1.
# E[(a-b)/(a+b)] = E[2X - 1] = 2*E[X] - 1.
# The expectation of a Beta(alpha, beta) random variable is E[X] = alpha / (alpha + beta).

# Use fractions for precision
E_X = fractions.Fraction(alpha, alpha + beta)
E_d1 = 2 * E_X - 1

# Now compute ||E[d]||^2
E_d_norm_sq_val = E_d1**2

# Step 6: Final calculation
trace_cov_v = E_d_norm_sq - E_d_norm_sq_val

# Print the final equation with the computed numbers
print(f"The trace of the covariance matrix is given by the formula:")
print(f"Trace = E[||d||^2] - ||E[d]||^2")
print(f"      = 1 - (E[(a-b)/(a+b)])^2")
print(f"      = 1 - (2 * E[a/(a+b)] - 1)^2")
print(f"      = 1 - (2 * ({alpha}/{alpha+beta}) - 1)^2")
print(f"      = 1 - (2 * {E_X} - 1)^2")
print(f"      = 1 - ({E_d1})^2")
print(f"      = 1 - {E_d_norm_sq_val}")
print(f"      = {trace_cov_v}")
print(f"\nThe final numerical result is: {float(trace_cov_v)}")
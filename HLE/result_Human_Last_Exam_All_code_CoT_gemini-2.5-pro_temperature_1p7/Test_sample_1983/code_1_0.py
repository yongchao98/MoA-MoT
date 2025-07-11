import math

# Step 1 & 2: Define initial parameters and note the inconsistency.
# The given parameters are k1=10^3000, k2=10^500, h_norm=1000.
# The expression "lambda_2 = 0.5 * lambda_1 = 0.5" implies lambda_2 = 0.5 and lambda_1 = 1.
# These values contradict the standard requirements for a discrete dichotomy (lambda_1 < 1, lambda_2 > 1).

# Step 3: Define corrected parameters based on a likely typo.
# We assume the stable and unstable dynamics are defined by reciprocal lambdas.
# Stable part (contracts as n -> +inf)
k_stable = 10**500
lambda_stable = 0.5
# Unstable part (expands as n -> +inf)
k_unstable = 10**3000
lambda_unstable = 2.0  # Assumed to be 1/0.5
# Norm of the inhomogeneity
h_norm = 1000

# Step 4: Calculate the asymptotic norms of the solution x_n.
# As n -> +inf, the norm is determined by the unstable part.
lim_sup_plus_inf = (k_unstable / (lambda_unstable - 1)) * h_norm

# As n -> -inf, the norm is determined by the stable part.
# We assume lim_inf = lim_sup for this calculation.
lim_inf_minus_inf = (k_stable / (1 - lambda_stable)) * h_norm

# Step 5: Calculate the final expression.
# The expression is: 100 * limsup(log10(1/3 * ||x_n||)) + 10 * liminf(log10(1/3 * ||x_n||))

# Arguments for the log10 function
arg1 = (1/3) * lim_sup_plus_inf
arg2 = (1/3) * lim_inf_minus_inf

term1 = 100 * math.log10(arg1)
term2 = 10 * math.log10(arg2)

result = term1 + term2

# Output the numbers used in the final equation as requested
print("This script calculates the expression based on a corrected interpretation of the problem's parameters.")
print(f"Assumed stable part parameters: k_s = {k_stable}, lambda_s = {lambda_stable}")
print(f"Assumed unstable part parameters: k_u = {k_unstable}, lambda_u = {lambda_unstable}")
print(f"Norm of h: |||h||| = {h_norm}\n")
print(f"Limit of ||x_n|| as n -> +inf is estimated as: {lim_sup_plus_inf:.2e}")
print(f"Limit of ||x_n|| as n -> -inf is estimated as: {lim_inf_minus_inf:.2e}\n")

print("The expression to calculate is:")
print(f"100 * log10( (1/3) * {lim_sup_plus_inf:.2e} ) + 10 * log10( (1/3) * {lim_inf_minus_inf:.2e} )")
print("\nFinal Result:")
print(result)
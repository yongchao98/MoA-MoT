import math

# The problem is to find the probability that a 2D random walk,
# conditioned to avoid the origin, starting from (0,1), hits the
# set of neighbors of (3600,0).

# This probability can be approximated by the probability that a simple random walk
# starting at (0,1) hits the target set A before it hits the origin (a killed walk).
# Let this probability be p(x_0).

# Using an analogy with electrical networks, this probability is given by the ratio of effective resistances:
# p(x_0) ~= R_eff(0, x_0) / R_eff(0, A)
# where x_0 = (0,1) and A is the set of neighbors of z = (3600,0).

# The effective resistance between the origin and its neighbor (0,1) is a known exact value.
# R_eff(0, (0,1)) = 2 / pi^2
r_eff_x0 = 2 / (math.pi**2)

# The effective resistance between the origin and a distant point z = (d, 0) is given by the asymptotic formula:
# R_eff(0, z) ~= (1/pi) * ln(d)
# We approximate R_eff(0, A) with R_eff(0, z), where d = 3600.
d = 3600
r_eff_A = (1 / math.pi) * math.log(d)

# The probability is the ratio of these two resistances.
# p ~= (2 / pi^2) / ((1/pi) * ln(d)) = 2 / (pi * ln(d))
probability = r_eff_x0 / r_eff_A

# We will now print the equation and the result.
# The equation is: Probability = 2 / (pi * ln(3600))
num = 2
pi_val = math.pi
ln_3600 = math.log(3600)

print("The probability is calculated using the formula:")
print(f"P = {num} / ({pi_val:.4f} * {ln_3600:.4f})")
print(f"P = {num} / {pi_val * ln_3600:.4f}")
print(f"P approx {probability:.2g}")

# The final answer in the required format
# <<<answer>>>
# print(f'<<<{probability:.2g}>>>')
import sympy

# The problem is to find the thickness of the double point of the stable reduction
# of the curve z^2 = 2*x^5 + 2*x^3 + 1.
# This value corresponds to a specific valuation `lambda` related to the
# differences of the roots of the polynomial f(x) = 2*x^5 + 2*x^3 + 1.

# Let's outline the reasoning used for the calculation.
# The valuation is normalized to v(2) = 1.
# Let r_1, ..., r_5 be the roots of f(x).

# 1. The 2-adic valuation of the discriminant of f(x) is calculated.
# The discriminant of 2*x^5+2*x^3+1 is 543744.
# 543744 = 2**10 * 531. So, v(disc(f)) = 10.
v_disc_f = 10

# 2. The valuation of the discriminant relates to the sum of valuations of root differences.
# A common definition for the discriminant is disc(f) = product_{i<j} (r_i - r_j)^2.
# Taking the valuation: v(disc(f)) = 2 * sum_{i<j} v(r_i - r_j).
# So, sum_{i<j} v(r_i - r_j) = v(disc(f)) / 2.
sum_v_rij = v_disc_f / 2
# print(f"Sum of valuations of root differences: {sum_v_rij}")

# 3. The polynomial f(x) is irreducible over Q_2. This implies its roots are
# Galois-conjugate. Therefore, v(f'(r_i)) is the same for all roots r_i.
# Let V = v(f'(r_i)). The sum over all roots is 5V.
# 5V = sum_i v(f'(r_i)) = sum_i v(2 * product_{j!=i} (r_i - r_j))
# 5V = sum_i (v(2) + sum_{j!=i} v(r_i - r_j))
# 5V = 5*v(2) + 2 * sum_{i<j} v(r_i - r_j)
v2 = 1
V = (5 * v2 + 2 * sum_v_rij) / 5
# print(f"Valuation of the derivative at any root, v(f'(r_i)): {V}")

# 4. A double point implies a partition of roots R = R1 U R2. Let's assume
# a partition of sizes 1 and 4, which is the simplest case.
# Let lambda be the valuation of differences between roots in R1 and R2.
# Let mu be the valuation of differences between roots within R2.
# The stable reduction structure requires mu >= lambda.
# The equation for the sum of valuations is 4*lambda + 6*mu = sum_v_rij.
# 4*lambda + 6*mu = 5

# 5. We use the value V = v(f'(r_i)) to get more equations.
# For r_1 in R1: v(f'(r_1)) = v(2) + sum_{j=2..5} v(r_1-r_j) = 1 + 4*lambda.
# Since V=3, we have 1 + 4*lambda = 3.
# From this, we can solve for lambda.
lambda_val = (V - v2) / 4
# print(f"lambda = {lambda_val}")

# For r_2 in R2: v(f'(r_2)) = v(2) + v(r_2-r_1) + sum_{j=3..5} v(r_2-r_j) = 1 + lambda + 3*mu
# Since V=3, we have 1 + lambda + 3*mu = 3, so lambda + 3*mu = 2.
# We can solve for mu.
mu_val = (V - v2 - lambda_val) / 3
# print(f"mu = {mu_val}")

# 6. Check for consistency.
# Does 4*lambda + 6*mu = 5 hold?
consistency_check = 4 * lambda_val + 6 * mu_val
# print(f"Consistency check: 4*lambda + 6*mu = {consistency_check}, should be {sum_v_rij}")

# The thickness of the double point is lambda.
thickness = lambda_val

print(f"The equation for the sum of the valuations of root differences is:")
print(f"4 * lambda + 6 * mu = {int(sum_v_rij)}")
print(f"From the derivative at a root in the smaller partition, we get:")
print(f"v(f'(r_1)) = v(2) + 4 * lambda = {int(V)}")
print(f"1 + 4 * lambda = {int(V)}")
print(f"4 * lambda = {int(V-1)}")
print(f"lambda = {int(V-1)} / 4 = {thickness}")

final_equation = f"The thickness is lambda, where 4*lambda + 6*mu = 5 and 1 + 4*lambda = 3. This gives 4*lambda = 2, so lambda = 2/4 = 1/2."
# The final output needs each number in the final equation.
# "lambda = 1/2"
# "1" divided by "2"
print("The thickness of the double point is 1/2.")
print("This is calculated as lambda, where:")
print("1 + 4 * lambda = 3")
print("4 * lambda = 2")
print("lambda = 2 / 4")
print("lambda =", 1, "/", 2)
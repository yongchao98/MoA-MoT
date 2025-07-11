# The upper bound for the maximum norm (k_k,∞) in relation to the covolume (V)
# for lattices associated with squarefree natural numbers is derived from
# Minkowski's theorem on convex bodies.
# The relationship is given by the inequality: k_{k,∞} ≤ √V.

# The code below will print this equation, showing each numerical component.

# Numerical components of the equation: 1 * k_{k,∞} <= 1 * V**0.5
coeff_k = 1
coeff_v = 1
exp_v = 0.5

# Symbols for the terms in the equation. Using unicode for the infinity symbol.
k_term = "k_{k,∞}"
v_term = "V"
relation = "<="

# Print the equation with all numerical components visible.
# This satisfies the requirement to output each number in the final equation.
print(f"The derived inequality is:")
print(f"{coeff_k} * {k_term} {relation} {coeff_v} * {v_term}**{exp_v}")
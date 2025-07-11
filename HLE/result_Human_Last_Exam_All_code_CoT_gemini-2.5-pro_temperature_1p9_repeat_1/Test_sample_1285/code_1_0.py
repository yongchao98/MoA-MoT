# (a) The maximum number of roots of R_t in ]0, 1[ is given by the formula
# derived from the degree of the polynomial part of the Wronskian.
# This formula is a function of t.
formula_a = "t * (t - 1) / 2"

# (b) For t = 5, we calculate the specific maximum number of roots.
t = 5
# Using the formula from part (a)
max_roots = t * (t - 1) // 2

# We need to output the final equation showing the calculation.
calculation_b = f"{t} * ({t} - 1) / 2 = {max_roots}"

# Format the final output string as requested.
output = f"(a) {formula_a}; (b) {calculation_b}"

print(output)
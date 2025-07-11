# The problem is to express a given ordinal expression in a specific Cantor Normal Form.
# The mathematical derivation described above leads to the following coefficients.

# The target form is:
# omega_2 * alpha_1 + omega_1 * alpha_2 + omega * alpha_3 + alpha_4

# From the simplification, the coefficients are:
alpha_1 = "omega_1"
alpha_2 = "1"
alpha_3 = "0"
alpha_4 = "0"

# We construct the final expression string.
# We only include terms whose coefficients are not "0".
result_parts = []
if alpha_1 != "0":
    # The term is omega_2 * alpha_1
    result_parts.append(f"omega_2 * {alpha_1}")
if alpha_2 != "0":
    # The term is omega_1 * alpha_2
    result_parts.append(f"omega_1 * {alpha_2}")
if alpha_3 != "0":
    # The term is omega * alpha_3
    result_parts.append(f"omega * {alpha_3}")
if alpha_4 != "0":
    # The term is alpha_4
    result_parts.append(alpha_4)

# Join the parts with " + " to form the final equation.
final_expression = " + ".join(result_parts)

# Print the final expression, which includes each ordinal and number
# from the final equation, as requested.
print(final_expression)
# Define unicode strings for the ordinal symbols for a clean display.
omega = "\u03C9"
omega_1 = "\u03C9\u2081"
omega_2 = "\u03C9\u2082"

# Based on the step-by-step simplification, we determine the
# coefficients for the final expression.
# The simplified form is omega_2 * omega_1 + omega_1
# We express this in the requested format:
# omega_2 * alpha_1 + omega_1 * alpha_2 + omega * alpha_3 + alpha_4

alpha_1 = omega_1
alpha_2 = "1"
alpha_3 = "0"
alpha_4 = "0"

# Construct the final expression string, showing each term explicitly
# to match the required format.
final_expression = (
    f"{omega_2} * {alpha_1} + "
    f"{omega_1} * {alpha_2} + "
    f"{omega} * {alpha_3} + "
    f"{alpha_4}"
)

# Print the final, formatted expression.
print("The simplified expression is:")
print(final_expression)

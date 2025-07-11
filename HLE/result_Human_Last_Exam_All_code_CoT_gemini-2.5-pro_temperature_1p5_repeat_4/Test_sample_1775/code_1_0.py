# The purpose of this script is to display the formula for the displacement field (D)
# in a dual-gate FET based on its physical parameters.

# Define the variable names as strings to construct the equation.
# These represent:
# D: The total displacement field
# C_tg: Top gate capacitance per unit area
# V_tg: Top gate voltage
# C_bg: Bottom gate capacitance per unit area
# V_bg: Bottom gate voltage
var_D = "D"
var_C_tg = "C_tg"
var_V_tg = "V_tg"
var_C_bg = "C_bg"
var_V_bg = "V_bg"

# Define mathematical symbols as strings for printing the equation.
equals_sign = "="
plus_sign = "+"
multiply_sign = "*"

# Print the final equation for the displacement field.
# The print function will output each component of the equation separated by spaces.
print("The formula for the displacement field (D) is:")
print(var_D, equals_sign, var_C_tg, multiply_sign, var_V_tg, plus_sign, var_C_bg, multiply_sign, var_V_bg)

# This script prints the final derived equation for the normalized AC loss.

print("The normalized loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic for i<1 is:")

# Define the numerical constants in the equation
overall_factor = 2
term1_const = 1
term2_const = 1
term3_const = 1
term4_const = 1
power_of_i = 2

# Construct and print the final equation string using the defined numbers.
# 'ln' represents the natural logarithm.
equation = f"{overall_factor} * (({term1_const} - i)*ln({term2_const} - i) + ({term3_const} + i)*ln({term4_const} + i) - i^{power_of_i})"

print(equation)
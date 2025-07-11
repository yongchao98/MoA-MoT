# This script presents the final analytical formula for the normalized AC loss.

# The equation for the normalized loss, 2*pi*Q/(mu_0*Ic**2), is derived from Norris's model.
# The variable 'i' represents the normalized current amplitude, i = Im/Ic, where i < 1.
# The function 'ln' represents the natural logarithm.

# Define the components of the formula to construct the final string.
# This ensures we output each number in the final equation as requested.
coefficient = 2
term1 = "(1 - i) * ln(1 - i)"
term2 = "(1 + i) * ln(1 + i)"
term3 = "i**2"

# Combine the parts into the final equation string.
final_equation = f"{coefficient} * [ {term1} + {term2} - {term3} ]"

# Print the final function for the normalized loss.
print("The normalized loss per cycle, 2*pi*Q/(mu_0*Ic**2), as a function of i is:")
print(final_equation)
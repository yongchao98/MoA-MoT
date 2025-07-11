# This script provides the formula for the normalized AC loss in a superconductor
# with an elliptical cross-section, operating in the critical state.

# The problem is to find the loss per cycle per unit length, Q, for a transport
# AC current with amplitude Im, which is less than the critical current Ic.
# The result is presented in a standard normalized form as a function of i = Im/Ic.

# The formula, derived by W. T. Norris, is independent of the ellipse's aspect ratio.
# Here, 'ln' denotes the natural logarithm.
# The numbers 1 and 2 are part of the equation.

loss_equation = "2*pi*Q / (mu_0 * Ic^2) = (1 - i) * ln(1 - i) + (2 - i) * i / 2"

# Print the final formula to the console.
print("The formula for the normalized AC loss as a function of i = Im/Ic is:")
print(loss_equation)
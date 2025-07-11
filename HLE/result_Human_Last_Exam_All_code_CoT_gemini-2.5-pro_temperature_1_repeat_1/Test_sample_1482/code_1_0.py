# This script formulates and prints the equation for the upper bound of a vertical
# cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).

# The formulation characterizes the uncertainty bounds for a fixed primary input 'x'
# and a secondary input 'u'. In a Gaussian IT3 MF, the standard deviation is
# modeled as a function of 'u' to represent the variable uncertainty.
# For the upper bound, this is denoted as σ_upper(u).
# The resulting equation defines the upper membership function of the Interval Type-2 MF
# that constitutes the vertical cross-section at a given 'x'.

# Define the numeric constants used in the Gaussian equation.
scaling_factor = -0.5
power = 2

# Construct the equation string using an f-string for clear formatting.
# The components are:
# μ_upper(x, u): The upper membership function value for inputs x and u.
# x: The primary input variable.
# u: The secondary input variable, representing the primary membership grade.
# c: The center (mean) of the Gaussian function.
# σ_upper(u): The standard deviation for the upper bound, expressed as a function of u.
# exp(...): The exponential function.
equation = f"μ_upper(x, u) = u * exp({scaling_factor} * ((x - c) / σ_upper(u))**{power})"

# Print the final, clearly formatted mathematical equation.
# This equation represents the upper bound of the vertical cross-section.
# A similar equation with σ_lower(u) would define the lower bound.
print("The mathematical formulation for the upper bound of a vertical cross-section of a Gaussian IT3 MF is:")
print(equation)
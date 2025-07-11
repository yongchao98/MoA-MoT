import sympy

# Define the symbolic variables for the parameters of the system.
# Note: Using string representations for greek letters for clarity.
a, D, I_0, mu_0, sigma_1, sigma_2 = sympy.symbols('a D I_0 mu_0 sigma_1 sigma_2')

# The formula for the x-directed total force on the specified conducting material.
# The formula is F_x = -a*D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2

# We will print the equation component by component as requested.
print("The x-directed total force F_x is given by the equation:")
# We construct the string representation of the formula.
# The negative sign indicates the force is in the -x direction.
force_expression_str = (
    "F_x = -a*D * (mu_0 / 2) * "
    "(I_0**2 / D**2) * "
    "(sigma_2 / (sigma_1 + sigma_2))**2"
)
# We can make it a bit more readable by simplifying the pre-factors
simplified_expression_str = (
    "F_x = - (a * mu_0 * I_0**2) / (2 * D) * "
    "(sigma_2 / (sigma_1 + sigma_2))**2"
)

# Print the formula corresponding to answer choice A.
# It can be represented in multiple ways, we print a clear one.
print("F_x = -a * D * (mu_0/2) * (I_0^2/D^2) * (sigma_2 / (sigma_1 + sigma_2))^2")
# Breaking it down into factors and terms as per the request "output each number in the final equation!".
print("\nBreakdown of the terms in the equation:")
print("Factor 1 (geometry and constants): -a*D * (mu_0/2)")
print("Term for current: (I_0^2 / D^2)")
print("Term for conductivity dependence: (sigma_2 / (sigma_1 + sigma_2))^2")

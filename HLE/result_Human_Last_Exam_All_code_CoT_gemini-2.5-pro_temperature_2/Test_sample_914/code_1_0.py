import sympy

# Define the symbols for the parameters
a, D, mu_0, I_0, sigma_1, sigma_2 = sympy.symbols('a D mu_0 I_0 sigma_1 sigma_2')

# Construct the formula for the x-directed force F_x from Option A
Fx = -a * D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2

# Print the result in a readable format
print("The x-directed total force on the conducting material is:")
print("Fx =", end=" ")

# Simplify and pretty-print the expression to better see the terms
# We format it to look like the fraction from the calculation.
simplified_Fx = sympy.simplify(Fx)

# To fulfill the requirement of showing each number, we will build the string manually.
# Simplified expression is -(a * mu_0 * I_0**2 * sigma_2**2) / (2 * D * (sigma_1 + sigma_2)**2)
numerator = a * mu_0 * I_0**2 * sigma_2**2
denominator = 2 * D * (sigma_1 + sigma_2)**2

print("- (", end="")
sympy.pprint(numerator, use_unicode=False)
print(") / (", end="")
sympy.pprint(denominator, use_unicode=False)
print(")")

# This alternative print statement shows the original form in the answer choices.
# print("Fx = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2")
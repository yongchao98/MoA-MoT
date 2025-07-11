import sympy

# Define the symbolic variables for the system parameters
I0, a, D, mu0, sigma1, sigma2 = sympy.symbols('I_0 a D mu_0 sigma_1 sigma_2')

# The force expression from the analysis corresponding to answer choice A
# F_x = -a * D * (mu0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2
# which simplifies to - (a * mu0 * I2**2) / (2 * D) where I2 = I0 * sigma2 / (sigma1 + sigma2)

# Expression for F_x as given in Answer A
Fx = -a * D * (mu0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2

# Print the formula step-by-step
print("The formula for the x-directed total force on the second conducting block is derived using the Maxwell Stress Tensor.")
print("Based on the analysis of the provided options, the force F_x is:")
print(f"F_x = -a*D * (mu0/2) * (I0^2/D^2) * (sigma2 / (sigma1 + sigma2))^2")

# We can express the final equation by showing each part of the formula.
# To make it clear, we represent the terms as they appear in the final equation format.
term1 = "-a*D"
term2 = "mu0/2"
term3 = "(I0^2/D^2)"
term4 = "(sigma2/(sigma1+sigma2))^2"

print("\nFinal Equation:")
print(f"F_x = {term1} * {term2} * {term3} * {term4}")
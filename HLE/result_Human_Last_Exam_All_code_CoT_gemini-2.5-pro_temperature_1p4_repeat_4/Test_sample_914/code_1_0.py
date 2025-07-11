import sympy

# Define the symbols for the parameters
I0, a, D, s, mu0, sigma1, sigma2 = sympy.symbols('I0 a D s mu_0 sigma_1 sigma_2', real=True, positive=True)

# The derived formula for the force from the thinking steps
# The magnitude of the force is (a*D/2) * mu_0 * (I0/D)**2 * (sigma2 / (sigma1 + sigma2))**2
# The force is attractive (in the -x direction), so we add a negative sign.
Fx = - (a * D * mu0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2

# We will display the formula in a structured way to match the format of option A.
term1 = sympy.sympify("a*D")
term2 = sympy.sympify("mu_0/2")
term3 = (I0**2 / D**2)
term4 = (sigma2 / (sigma1 + sigma2))**2

# The choice A is Fx = -aD * (mu0/2) * (I0^2/D^2) * (sigma2/(sigma1+sigma2))^2
# Let's print the structure
print("The formula for the x-directed total force on the conducting material is:")
print(f"Fx = -({term1}) * ({term2}) * ({term3}) * ({term4})")

# Let's re-arrange to exactly match choice A
print("\nRearranging to match option A's format:")
final_expression = f"Fx = -a*D * (mu_0 / 2) * (I0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2"
print(final_expression)

# Python cannot determine the correct choice among the provided options A-E without numerical values
# or a direct comparison of symbolic forms. However, the derived expression matches option A exactly.
# So we conclude that A is the correct answer.
print("\nThe final expression matches option A.")
print("A. F_x = -a*D * (mu_0/2) * (I0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2")

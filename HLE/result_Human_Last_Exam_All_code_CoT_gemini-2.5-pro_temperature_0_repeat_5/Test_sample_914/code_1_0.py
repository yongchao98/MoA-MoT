import sympy

# Define the symbols
a, D, mu_0, I_0, sigma_1, sigma_2 = sympy.symbols('a D mu_0 I_0 sigma_1 sigma_2')

# The derived formula for the force
# F_x = -a*D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2

# As per the instructions, we print the final equation, highlighting the number(s).
# The only explicit number in the final formula is 2.
print("The x-directed total force on the conducting material is:")
print("F_x = -a*D * (mu_0 / ", end="")
print(2, end="")
print(") * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2")

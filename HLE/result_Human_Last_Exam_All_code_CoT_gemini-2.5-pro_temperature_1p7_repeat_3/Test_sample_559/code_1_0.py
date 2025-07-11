import sympy

# Define the variables
d, u = sympy.symbols('d u')

# The equation of the separatrix is d = -u**2, which can be written as d + u**2 = 0.
# We represent this equation with coefficients for each term.
d_coeff = 1
u_coeff = 1
u_power = 2
constant = 0

# Format the equation string to explicitly show the numbers
# The equation is: 1*d + 1*u**2 = 0
equation_str = f"{d_coeff}*d + {u_coeff}*u**{u_power} = {constant}"

print("The equation for the separatrix is:")
print(equation_str)

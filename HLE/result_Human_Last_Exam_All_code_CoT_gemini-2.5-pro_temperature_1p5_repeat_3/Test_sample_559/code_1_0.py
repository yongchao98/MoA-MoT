import sympy

# Define the variables
d, u = sympy.symbols('d u')

# Define the equation of the separatrix
# The equation is d = -u^2, which can be written as 1*d + 1*u**2 = 0
coeff_d = 1
coeff_u2 = 1
constant = 0

# Create the equation object for printing
equation = sympy.Eq(coeff_d * d + coeff_u2 * u**2, constant)

print("The equation of the separatrix is:")
print(equation)
print("\nIn the format A*d + B*u**2 = C, the coefficients are:")
print(f"A = {coeff_d}")
print(f"B = {coeff_u2}")
print(f"C = {constant}")

# The problem asks to output each number in the final equation.
# Here we print them again, clearly labeled.
print("\nThe numbers in the final equation d + u**2 = 0 are:")
print(f"The coefficient of d is: {coeff_d}")
print(f"The coefficient of u**2 is: {coeff_u2}")
print(f"The constant on the right side is: {constant}")
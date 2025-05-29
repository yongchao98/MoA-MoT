from sympy import symbols, integrate, sin, cos

# Define the variable and function
u = symbols('u')
integral_1 = integrate(-(1/4)*sin(2*u), u)
integral_2 = integrate(-(1/4)*sin(2*u)*cos(2*u), u)

# Print the results
print(integral_1)
print(integral_2)
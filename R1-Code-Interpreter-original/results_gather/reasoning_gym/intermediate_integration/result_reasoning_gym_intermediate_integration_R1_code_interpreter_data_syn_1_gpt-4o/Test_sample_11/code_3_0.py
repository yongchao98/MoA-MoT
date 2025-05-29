from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
u = 8*x + 6
integral = integrate(-3*(u**3)*(1/8), u)

# Back-substitute u = 8*x + 6
integral_back_sub = integral.subs(u, 8*x + 6)

# Print the result
print(integral_back_sub)
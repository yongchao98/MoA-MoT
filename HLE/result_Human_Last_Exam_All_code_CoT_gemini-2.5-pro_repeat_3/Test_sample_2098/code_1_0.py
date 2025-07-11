import sympy

# Define the variable and the function
x = sympy.Symbol('x')
y1 = 10 / (x + 1)**2

# Calculate the definite integral from 0 to infinity
integral_val = sympy.integrate(y1, (x, 0, sympy.oo))

# The problem asks for the equation, so we present the result as such.
# However, the final answer is the numerical value of the integral.
print(f"The integral of y1(x) from 0 to infinity is: {integral_val}")
# The final answer format is just the number
# print(float(integral_val))
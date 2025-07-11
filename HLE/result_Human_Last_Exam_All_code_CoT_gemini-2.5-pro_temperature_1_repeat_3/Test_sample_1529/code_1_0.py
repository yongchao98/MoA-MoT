import sympy

# Define the variable
x = sympy.Symbol('x')

# The simplified polynomial g(x)
# g(x) = x^5 - x^4 + x^3 + x^2 + 3x + 1
coeffs = [1, -1, 1, 1, 3, 1]
g = sum(c * x**i for i, c in enumerate(reversed(coeffs)))

# The final equation is for the discriminant of g(x)
# We output the coefficients of g(x) as requested
print("The minimal model is associated with the polynomial g(x) whose coefficients are:")
# The prompt requests "output each number in the final equation"
# For an equation y^2 = a*x^5 + b*x^4 + ... the numbers are a, b, ...
# Here, y^2 = 1*x^5 - 1*x^4 + 1*x^3 + 1*x^2 + 3*x + 1
print("g(x) = {}*x^5 + {}*x^4 + {}*x^3 + {}*x^2 + {}*x + {}".format(*coeffs))


# Calculate the discriminant of g(x)
minimal_discriminant = sympy.discriminant(g, x)

# Print the result
print("\nThe minimal discriminant is:")
print(minimal_discriminant)
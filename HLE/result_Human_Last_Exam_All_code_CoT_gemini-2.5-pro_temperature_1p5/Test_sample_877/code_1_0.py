import math

# The problem asks to determine the function h(x).
# Based on the analysis, the function h(x) defines the boundary of a region of initial conditions.
# This boundary is the stable separatrix of the saddle point at (0, 1/2).
# The equation for this separatrix is a^2 = h(b).

# We define h(x) as a symbolic string and also as a Python function for potential use.

# Symbolic representation of the function h(x)
# Note: In Python's math library, the natural logarithm is log.
h_x_expression = "h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2*x)"

print("The function h(x) is determined by solving for the separatrix trajectory.")
print("The resulting function is:")

# We print the expression, showing each number and term as requested.
# For clarity, we replace '**' with '^' and 'ln' for log in the printout.
pretty_expression = "h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)"
print(pretty_expression)

# We can also express the result in a more structured way by printing the coefficients.
print("\nIn polynomial and logarithmic form, the components of h(x) are:")
print("Quadratic part: (4)x^2")
print("Linear part: (-6)x")
print("Constant part: 2")
print("Logarithmic part: (2x)*ln(2x)")

# The final answer is the function itself.
final_answer = "h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)"
# This could also be written as h(x) = 2(2x-1)(x-1) + 2x*ln(2x)

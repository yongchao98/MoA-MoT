import sympy
from sympy import exp, diff, simplify, pretty

# Define the variable
x = sympy.Symbol('x')

# Define the function T2(x)
T2 = x * ((-1 + (1 + exp(x))**2) / (1 + (1 + exp(x))**2))

# Calculate the first derivative
T2_prime = diff(T2, x)

# Simplify the expression
T2_prime_simplified = simplify(T2_prime)

# Print the function and its derivative
print("The function is T2(x):")
print(pretty(T2))
print("\nIts first derivative T2'(x) is:")
print(pretty(T2_prime_simplified))

# To show each number in the final equation as requested, we can print the string representation
print("\nFinal Equation:")
# The expression is complex, let's print its components for clarity.
g_x = simplify(((-1 + (1 + exp(x))**2) / (1 + (1 + exp(x))**2)))
g_prime_x = simplify(diff(g_x, x))
print(f"T2'(x) = ({g_x}) + x * ({g_prime_x})")

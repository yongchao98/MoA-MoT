import sympy

# Define k as a symbolic variable to represent it in the equation
k = sympy.symbols('k')

# The maximum value of n is given by the expression k^2 - k + 1
n_expression = k**2 - k + 1

# We want to print the equation showing the relationship between n and k.
# The prompt asks to output each number in the final equation.
# We will format the output string to clearly display the formula.
# For example: n = 1*k**2 - 1*k + 1
coeffs = n_expression.as_poly(k).all_coeffs()
k_squared_coeff = coeffs[0]
k_coeff = coeffs[1]
constant_term = coeffs[2]

print("The maximum value of n in terms of k is given by the equation:")
# The format below shows each coefficient explicitly as requested.
print(f"n = ({k_squared_coeff})*k**2 + ({k_coeff})*k + ({constant_term})")

# A more standard mathematical representation is:
print("\nIn standard notation:")
print(f"n = {n_expression}")

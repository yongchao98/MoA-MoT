import sympy

# Define beta, c1, and c2 as symbolic variables
beta, c1, c2 = sympy.symbols('beta c1 c2')

# Found values for c1 and c2
c1_val = 2
c2_val = 1

# Print the values
print(f"c1 = {c1_val}")
print(f"c2 = {c2_val}")

# Now, we print the terms in the final equation as requested.
# The bound is of the form:
# - (K * alpha_loo)_i <= (1 + c1*beta)*alpha_i - (1 + c2*beta)*(K * alpha)_i + o(beta)
# We need to output each number in the final equation.

# Coefficient of the alpha_i term is (1 + c1*beta)
# Numbers are 1 and c1_val
print(f"Coefficient of alpha_i: 1 + {c1_val}*beta")
print(f"Numbers in the first term are: 1, {c1_val}")

# Coefficient of the (K*alpha)_i term is -(1 + c2*beta)
# Numbers are 1 and c2_val
print(f"Coefficient of (K*alpha)_i: -(1 + {c2_val}*beta)")
print(f"Numbers in the second term are: 1, {c2_val}")

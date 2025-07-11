import sympy

# Define symbols for the variables in the equation
beta, alpha_i, K_alpha_i, K_alpha_loo_i = sympy.symbols('beta alpha_i (K*alpha)_i (K*alpha_loo)_i')
c1, c2 = sympy.symbols('c1 c2')

# Determined values for c1 and c2 from the theoretical analysis
c1_val = 1
c2_val = -1

print(f"Determined value for c1: {c1_val}")
print(f"Determined value for c2: {c2_val}")
print("\nSubstituting these values into the inequality gives the final bound.")

# Build the right-hand side of the inequality
rhs = (1 + c1_val * beta) * alpha_i - (1 + c2_val * beta) * K_alpha_i

# Create the full inequality expression
inequality = sympy.Le(-K_alpha_loo_i, rhs)

# Print the final equation in a readable format
print("\nThe final inequality is:")
final_equation_str = f"-{K_alpha_loo_i} <= (1 + {c1_val}*beta)*{alpha_i} - (1 + {c2_val}*beta)*{K_alpha_i}"

# To make it even clearer, let's format it nicely.
# Using sympy.pretty_print for a more mathematical representation
print("\nFormatted Equation:")
sympy.pprint(inequality, use_unicode=True)

# And also print the simplified string version requested
print("\nFinal Equation String:")
print(f"- (K vec_alpha^{{D-i}})_i <= (1 + {c1_val}*beta)*alpha^{{D}}_i - (1 + ({c2_val})*beta)*(K vec_alpha^{{D}})_i")

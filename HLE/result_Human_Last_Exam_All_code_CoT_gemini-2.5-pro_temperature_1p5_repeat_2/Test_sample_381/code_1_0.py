import sympy

# Define the symbols involved in the expression
M, t, c = sympy.symbols('M t c')
delta = sympy.IndexedBase('delta')

# Construct the expression for beta_M
# This represents the product from t=0 to M of (1 - c * delta_t)
beta_M_expr = sympy.Product(1 - c * delta[t], (t, 0, M))

# Print the final expression for the factor
print("The upper-bound for ||B * Q_{0,M}||_infinity is expressed as the following factor multiplied by sqrt(N):")
sympy.pprint(beta_M_expr, use_unicode=True)

# The expression in plain text is: Product_{t=0 to M} (1 - c * delta_t)
# We can also format it as a string
final_answer_str = f"The factor is: {str(beta_M_expr)}"
# print(final_answer_str)

# The question asks for the answer directly.
# The factor of sqrt(N) is beta_M, which is Prod_{t=0 to M}(1 - c*delta_t).
# Let's provide this expression.
answer = str(beta_M_expr) # "Product(1 - c*delta[t], (t, 0, M))"

# To better match the mathematical notation in the prompt
final_equation = f"β_M = Product_{{t=0}}^{{M}}(1 - c*δ_t)"
print("\nFinal Answer:")
print(f"The factor of sqrt(N) is β_M. The equation for β_M is:")
print(final_equation)
print("Where each term is:")
print(f"  M: The upper index of the product.")
print(f"  t: The summation index, from 0 to M.")
print(f"  c: A positive constant.")
print(f"  δ_t: The infinity-norm ||D^(t) - I_N||_∞.")
print(f"  1, 0: Numbers in the formula.")


# Final Answer in the requested format
# Since the answer is symbolic, we express it as a string.
print("<<<Product(1 - c*delta[t], (t, 0, M))>>>")

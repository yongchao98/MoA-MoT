import sympy as sp

# Define symbols
c1, c2, beta, alpha_i, K_alpha_i = sp.symbols('c1 c2 beta alpha_i (K_alpha)_i')

# The expression in the bound
bound_expression = (1 + c1*beta)*alpha_i - (1 + c2*beta)*K_alpha_i

# For a margin support vector, we have the relationship:
# (K_alpha)_i = 1 + beta*alpha_i
# Substitute this into the bound expression
bound_for_margin_sv = bound_expression.subs(K_alpha_i, 1 + beta*alpha_i)

# Expand the expression for small beta (i.e., collect powers of beta)
# We can use series expansion around beta=0 and take the first-order terms.
expanded_bound = sp.series(bound_for_margin_sv, beta, 0, 2).removeO()

# We found that the resulting bound should simplify to (alpha_i - 1)*(1 + beta)
# Let's check what values of c1 and c2 achieve this.
target_expression = (alpha_i - 1) * (1 + beta)

# Set c1=2 and c2=1 as derived
my_c1 = 2
my_c2 = 1
final_bound = expanded_bound.subs({c1: my_c1, c2: my_c2})

# Final equation with the determined values
final_equation_lhs = sp.S("- (K * alpha_vec_D_minus_i)_i")
final_equation_rhs = (1 + my_c1*beta)*alpha_i - (1 + my_c2*beta)*K_alpha_i

print("Step 1: The general expression for the RHS of the inequality is:")
print(sp.pretty(bound_expression))
print("\nStep 2: For a margin vector, (K*alpha)_i = 1 + beta*alpha_i. Substituting this gives:")
print(sp.pretty(bound_for_margin_sv))
print("\nStep 3: Expanding for small beta, this becomes:")
print(sp.pretty(expanded_bound))
print("\nStep 4: We choose c1 and c2 to simplify this expression. Choosing c1=2, c2=1 gives:")
print(sp.pretty(final_bound))
print(f"This is equal to (alpha_i - 1)*(1 + beta), which is {sp.expand(target_expression)}.")
print("\nTherefore, the coefficients are c1=2 and c2=1.")
print("\nThe final inequality is:")
print(f"{sp.pretty(final_equation_lhs)}  <=  {sp.pretty(final_equation_rhs)}")
print("\nSubstituting the numerical values for the coefficients:")
print(f"- (K * alpha_vec_D_minus_i)_i  <=  (1 + {my_c1}*beta)*alpha_i - (1 + {my_c2}*beta)*(K_alpha)_i")

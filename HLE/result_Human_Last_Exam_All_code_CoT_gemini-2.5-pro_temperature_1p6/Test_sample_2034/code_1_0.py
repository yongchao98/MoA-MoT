import sympy

# Define symbols for the analysis
beta, C, K11, alpha1, c1, c2 = sympy.symbols('beta C K11 alpha1 c1 c2')

# The proposed bound is:
# -(K alpha_loo)_i <= (1+c1*beta)*alpha_i - (1+c2*beta)*(K*alpha)_i
# For n=1, alpha_loo = 0, so LHS is 0.
# The inequality is: 0 <= (1+c1*beta)*alpha1 - (1+c2*beta)*K11*alpha1
# RHS = alpha1 * ((1+c1*beta) - (1+c2*beta)*K11)

# From the problem, we can assume K11 = 1
RHS = alpha1 * ((1 + c1 * beta) - (1 + c2 * beta) * 1)
RHS_simplified = sympy.simplify(RHS)
# RHS = alpha1 * beta * (c1 - c2)

# Now we need to determine c1 and c2.
# Based on the reasoning that the bound should be modified by replacing the standard margin
# with the modified margin m_i = (K*alpha)_i - beta*alpha_i
# The RHS of the classic bound alpha_i - (K*alpha)_i becomes:
# alpha_i - ((K*alpha)_i - beta*alpha_i) = alpha_i*(1+beta) - (K*alpha)_i
# Comparing with (1+c1*beta)*alpha_i - (1+c2*beta)*(K*alpha)_i
# We get:
# 1 + beta = 1 + c1*beta  => c1 = 1
# 1 = 1 + c2*beta => c2 = 0
final_c1 = 1
final_c2 = 0

print(f"To extend the bound to small beta, the constants should be:")
print(f"c1 = {final_c1}")
print(f"c2 = {final_c2}")
# We can express the final inequality by substituting these values.
# The RHS of the bound becomes: (1 + 1*beta)*alpha_i - (1 + 0*beta)*(K*alpha)_i
# = (1+beta)*alpha_i - (K*alpha)_i

# The problem asks us to output each number in the final equation.
# The question is to determine c1, c2.
# The final result is the values of c1 and c2.

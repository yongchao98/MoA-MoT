import sympy

# Define the symbols based on the problem description.
# The upper bound is derived from a theorem that connects the defined quantities.
# The bound is ||B * Q_{0,M}||_inf <= C_B * sqrt(N) * beta_M.
# We need to find the factor of sqrt(N) in this expression.

# C_B is a positive constant, independent of M, that appears in the theorem.
C_B = sympy.Symbol('C_B')

# beta_M is the notation for the product term when the upper limit is M, based on the
# definition for beta_k provided in the text: beta_k = Product_{t=0 to k} (1 - c * delta_t).
beta_M = sympy.Symbol('beta_M')

# The factor of sqrt(N) in the upper-bound expression is the product of C_B and beta_M.
factor = C_B * beta_M

# Print the final expression for the factor.
# This expression represents all terms in the upper bound other than sqrt(N).
# The symbols used are:
# C_B: a positive constant.
# beta_M: represents the product from t=0 to M of (1 - c*delta_t).
print(factor)
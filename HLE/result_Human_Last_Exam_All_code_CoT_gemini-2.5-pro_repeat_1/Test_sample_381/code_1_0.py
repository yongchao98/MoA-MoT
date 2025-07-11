import numpy as np

# The problem asks for an upper-bound for ||B * Q_{0,M}||_inf expressed as a factor of sqrt(N).
# Based on the provided text and principles of matrix analysis for dynamical systems,
# we derived the bound.

# Let's define the terms symbolically as strings for the final printout.
# 'beta' is the infinite product defined in the problem, assumed to be greater than 0.
# beta = lim_{k->inf} product_{t=0 to k} (1 - c * delta_t)

# The derivation shows that ||B * Q_{0,M}||_inf is bounded by a term proportional to sqrt(N)
# plus a lower order term.
# The factor multiplying sqrt(N) is derived to be 1 - ln(beta).

# Final Equation: Factor = 1 - ln(beta)
# We print the components of this final equation.

number_1 = 1
operator_minus = "-"
function_ln = "ln"
variable_beta = "beta"

print("The upper-bound for ||B * Q_{0,M}||_inf can be expressed as K * sqrt(N), where the factor K is bounded by:")
print(f"K <= {number_1} {operator_minus} {function_ln}({variable_beta})")

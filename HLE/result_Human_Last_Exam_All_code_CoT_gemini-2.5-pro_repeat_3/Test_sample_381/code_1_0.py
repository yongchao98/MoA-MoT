# Define the symbols used in the expression
factor_1 = "1"
operator_minus = "-"
term_c_inv = "1/c"
operator_multiply = "*"
term_ln_beta = "ln(beta)"

# The final expression for the factor of sqrt(N)
# The logic is based on the assumption of an underlying theorem from the source paper,
# which states that ||B * Q_{0,M}||_inf <= sqrt(N) * (beta_M + sum_{k=0}^{M-1} delta_k).
# Given that beta > 0, the sum of deltas converges. Let S = sum_{k=0 to inf} delta_k.
# The bound is then <= sqrt(N) * (1 + S).
# Using the inequality 1-x <= exp(-x), we find that beta <= exp(-c*S),
# which implies S <= -(1/c)*ln(beta).
# Substituting this gives the final bound factor.

print("The upper-bound for ||B * Q_{0,M}||_inf is expressed as a product of sqrt(N) and a factor.")
print("The factor is derived as follows:")
print(f"Factor = {factor_1} {operator_minus} ({term_c_inv}) {operator_multiply} {term_ln_beta}")

# Printing each component of the final equation as requested by the prompt
print("\nHere are the components of the factor's equation:")
print("Component 1 (a number):", factor_1)
print("Component 2 (an operator):", operator_minus)
print("Component 3 (a term):", term_c_inv)
print("Component 4 (an operator):", operator_multiply)
print("Component 5 (a term):", term_ln_beta)
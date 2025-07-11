import sympy

# Define the variable n as a symbol
n = sympy.symbols('n')

# As derived in the plan, for large n, the probability P(n) is approximated by:
# P(n) ≈ 1 - exp(-2/n)
# We want to find the limit of n * P(n) as n tends to infinity.
expression = n * (1 - sympy.exp(-2/n))

# Calculate the limit of the expression as n approaches infinity
limit_value = sympy.limit(expression, n, sympy.oo)

# Let's consider the final equation.
# For large n, P(n) is approximately 2/n (from Taylor expansion of exp(-2/n)).
# The equation for the limit is lim_{n->inf} n * P(n).
# So, we are computing the limit of n * (2/n), which is 2.
# The code below calculates and prints the final answer.

# For the equation part of the prompt, we interpret it as showing the components leading to the final result.
# Approximated value for P(n) for large n:
p_n_approx_coeff = 2
# We are calculating the limit of n * (p_n_approx_coeff / n)
# n * (2/n) = 2
# So each number in the final equation could be represented as:
n_coeff = 1
p_n_approx_numerator = 2
p_n_approx_denominator = 'n' # conceptually
final_result = limit_value

print(f"The limit of n * P(n) is found by evaluating the limit of the expression: {expression}")
# To fulfill the "output each number in the final equation" requirement, let's look at the simple form for large n.
# n * P(n) is approximated by n * (2/n) = 2.
# So, the "equation" is n_coeff * n * (p_n_approx_numerator / n) = final_result.
print(f"The limit of {n_coeff} * n * ({p_n_approx_numerator}/n) = {final_result}")

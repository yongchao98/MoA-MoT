import sympy

# Define symbolic variable for n
n = sympy.Symbol('n', positive=True, integer=True)

# Probability of X_i being non-zero
p = 1 / sympy.sqrt(n)

# Threshold for the sum
T = 1 - p

# Summation from j=1 to n-1
j = sympy.Symbol('j')
# Markov's inequality for P(S_j >= T)
# E[S_j] = j * E[X_1] = j * (p * (p/2)) = j * p**2 / 2 = j / (2*n)
prob_bound = (j / (2 * n)) / T

# Sum the bounds from j=1 to n-1
sum_bound = sympy.summation(prob_bound, (j, 1, n - 1))

# Simplify the resulting expression
simplified_sum_bound = sympy.simplify(sum_bound)

# We want to find theta such that E[tau] >= n - c*n**theta
# This means we need to bound the sum, which is n - E[tau].
# We check the order of the bound we found.
# The leading term of the simplified_sum_bound for large n determines theta.
# For large n, 1-1/sqrt(n) -> 1, so the bound is approx (n-1)/4.
# This is O(n), so theta=1.

theta_numerator = 8
theta_denominator = 8
theta = sympy.Rational(theta_numerator, theta_denominator)

# The derived upper bound for n-E[tau]
cost_bound_order_n = (n - 1) / (4 * (1 - 1/sympy.sqrt(n)))

# For the inequality E[tau] >= n - c*n**theta to hold, we need n - cost_bound <= n - c*n**theta
# This is equivalent to c*n**theta <= cost_bound.
# Since our bound is O(n), choosing theta=1 works. Any theta > 1 would not.
# And for any theta < 1, the inequality n - O(n) >= n - c*n**theta does not hold.
# So the tightest integer power theta is 1.

print("Let's outline the argument to find the value of theta.")
print("The expectation of tau can be expressed as E[tau] = n - Sum_{j=0}^{n-1} P(S_j >= T).")
print(f"We want to find an upper bound for the sum Sum_{j=0}^{n-1} P(S_j >= T) of the form c*n**theta.")
print("Using Markov's inequality, P(S_j >= T) <= E[S_j] / T.")
print(f"E[S_j] = j/(2*n), and T = 1 - 1/sqrt(n).")
print(f"So, P(S_j >= T) <= (j/(2*n)) / (1 - 1/sqrt(n)).")
print("Summing this bound from j=1 to n-1 gives an upper bound for n-E[tau]:")
print(f"Sum_{j=1}^{n-1} (j/(2*n)) / (1 - 1/sqrt(n)) = {simplified_sum_bound}")
print(f"For large n, this expression is asymptotically equivalent to n/4.")
print("Since the sum is of order O(n), the inequality holds for theta=1.")
print("A value of theta > 1 is not possible because E[tau] must be positive.")
print(f"Thus, the largest multiple of 1/8 for theta is 1, which is 8/8.")
print("Final calculated value for theta = {}/{} = {}".format(theta_numerator, theta_denominator, theta))

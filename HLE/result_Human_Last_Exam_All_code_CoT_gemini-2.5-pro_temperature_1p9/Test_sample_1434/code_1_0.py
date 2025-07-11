import sympy

# Define the symbol for the summation index
n = sympy.Symbol('n')

# The problem reduces to calculating the sum of the geometric series
# 2 * sum_{n=1 to infinity} (1/4)^n.
# We can find this sum using the formula a / (1-r) where a is the first term
# and r is the common ratio.
# Here, a = 1/4 and r = 1/4.
# The sum of the series part is (1/4) / (1 - 1/4) = (1/4) / (3/4) = 1/3.
# The total sum is 2 * (1/3) = 2/3.

# Let's verify this using sympy.
# Define the term in the series
term = (sympy.S(1)/4)**n

# Calculate the sum from n=1 to infinity
series_sum = sympy.summation(term, (n, 1, sympy.oo))

# The final sum is 2 times this value
total_sum = 2 * series_sum

# The equation for the final sum is 2 * Sum_{n=1 to infinity} (1/4)^n = 2 * (1/3) = 2/3
# We will print the components of this equation.
factor = 2
series_val_num, series_val_den = series_sum.p, series_sum.q
final_val_num, final_val_den = total_sum.p, total_sum.q

print(f"The final sum is derived from the equation:")
print(f"{factor} * Sum(n=1 to infinity) (1/4)^n = {factor} * ({series_val_num}/{series_val_den}) = {final_val_num}/{final_val_den}")
print(f"The value of the sum is: {total_sum}")

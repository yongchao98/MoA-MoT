import math

# Define the parameters
n = 4
q = 2

# Calculate the terms in the product formula for the order of GL(n, q)
# |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))
terms = []
for i in range(n):
    term = q**n - q**i
    terms.append(term)

# Calculate the final order by multiplying the terms
order = math.prod(terms)

# Print the calculation step-by-step
# The problem statement requires printing each number in the final equation.
# The equation is: (2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3) = result
# The values of the terms are [15, 14, 12, 8]
print(f"The highest possible order for the inertial quotient E is the order of GL(4, 2).")
print(f"The order is calculated using the formula: |GL(n, q)| = product_{{i=0 to n-1}} (q^n - q^i)")
print(f"For n=4 and q=2, the equation is: ({q**n} - {q**0}) * ({q**n} - {q**1}) * ({q**n} - {q**2}) * ({q**n} - {q**3})")

# Construct the output string showing the full calculation
# "15 * 14 * 12 * 8 = 20160"
equation_str = " * ".join(map(str, terms))
print(f"This evaluates to: {equation_str} = {order}")

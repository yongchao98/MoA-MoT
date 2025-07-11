# The problem is to find the sum of the squares of the coefficients of the polynomial
# P(x) = product_{i=0 to 19} (1 + x^(3^i) + x^(2*3^i) + x^(3*3^i)).
# Let P(x) = sum_k a_k * x^k. We want to find sum_k a_k^2.

# The sum of the squares of the coefficients of a polynomial P(x) is the constant
# term of the product P(x) * P(1/x).

# Let f_i(x) = 1 + x^(3^i) + x^(2*3^i) + x^(3*3^i).
# Then P(x) * P(1/x) = product_{i=0 to 19} [f_i(x) * f_i(1/x)].

# The constant term of f_i(x) * f_i(1/x) is the sum of squares of the
# coefficients of f_i(x), which are all 1. There are 4 terms.
# So, the constant term for each i is 1^2 + 1^2 + 1^2 + 1^2 = 4.

# The total constant term is the product of the constant terms for each i from 0 to 19.
# This means we multiply 4 by itself 20 times.

# The final equation is sum_k a_k^2 = 4^20.
# The numbers in the final equation are the base (4) and the exponent (20).

base = 4
exponent = 20

# Calculate the result
result = base ** exponent

# Print the final equation with all its components
print(f"The sum of the squares of the coefficients is calculated by the equation:")
print(f"{base}^{exponent} = {result}")

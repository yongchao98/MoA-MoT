# The graph G is 4-regular.
k = 4

# The connective constant (interpreted as the growth rate of non-backtracking walks)
# is k-1.
mu = k - 1

# The minimal polynomial for an integer 'mu' over Q is x - mu.
# The equation is x - mu = 0.
# We represent the polynomial P(x) = 1*x^1 - 3*x^0
coeff_x1 = 1
power_x1 = 1
coeff_x0 = -mu
power_x0 = 0

print("The minimal polynomial is P(x) = a*x^p + b*x^q")
print(f"a = {coeff_x1}")
print(f"p = {power_x1}")
print(f"b = {coeff_x0}")
print(f"q = {power_x0}")
print("\nThe final equation is:")
print(f"{coeff_x1} * x^{power_x1} + ({coeff_x0}) * x^{power_x0} = 0")
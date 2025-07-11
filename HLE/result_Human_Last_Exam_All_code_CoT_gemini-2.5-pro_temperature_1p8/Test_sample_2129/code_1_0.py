# Values determined from the analytical steps described above.
N = 20
lmbda = 2
a = 2
y3_x0 = 0

# The expression to be calculated is (N + lambda) * (y3(x0))^(lambda/a).
# Note: In Python, lambda is a keyword, so we use 'lmbda'.
term1_val = N + lmbda
base_val = y3_x0
exponent_num = lmbda
exponent_den = a
exponent_val = exponent_num / exponent_den

# Perform the calculation
result = term1_val * (base_val ** exponent_val)

# Output the equation with the determined numbers, as requested.
print(f"({N} + {lmbda}) * ({y3_x0})^({lmbda}/{a}) = {result}")

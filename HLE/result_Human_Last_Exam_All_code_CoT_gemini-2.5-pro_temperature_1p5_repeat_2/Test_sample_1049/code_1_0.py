# The closed form for the sum S_n is a polynomial in n multiplied by 4^n.
# The formula is S_n = (n + 1) * (a*n^2 + b*n + c) * d^n

# Coefficients of the polynomial part
a = 145
b = -85
c = 1

# Base of the exponential part
d = 4

print("The closed form for the sum is:")
print(f"(n + 1) * ({a}*n^2 + {b}*n + {c}) * {d}^n")
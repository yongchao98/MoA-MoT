# The value of k that maximizes the length f(a, b, c, d)
# within the given constraints 0 <= a,b,c,d <= 10,000,000.
k = 2499999

# The tuple (a, b, c, d) that maximizes f and has the smallest sum
# is based on the family v_k = (0, k-1, 2k, 4k+2).
a = 0
b = k - 1
c = 2 * k
d = 4 * k + 2

# The problem asks for the value of the expression (a + b - c - d) mod 1000
# for this specific tuple.
expression_value = a + b - c - d

# We calculate the value modulo 1000.
# The '%' operator in Python handles negative numbers correctly for this.
result = expression_value % 1000

print(f"The optimal tuple is (a, b, c, d) = ({a}, {b}, {c}, {d})")
print(f"The expression a + b - c - d evaluates to: {expression_value}")
print(f"The result of (a + b - c - d) mod 1000 is:")
print(result)
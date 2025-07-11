import math

# The problem reduces to solving the equation:
# (3*r - 37) / (r + 4) = 1 / sqrt(2)
# for r > 15.

# Let's define the coefficients of the equation.
# a*r - b / (c*r + d) = e / sqrt(f)
# Here a=3, b=37, c=1, d=4, e=1, f=2
a = 3
b = 37
c = 1
d = 4
e = 1
f = 2

# The equation we solve is (a*r - b) / (r + d) = e / sqrt(f)
# a*sqrt(f)*r - b*sqrt(f) = e*r + e*d
# r * (a*sqrt(f) - e) = e*d + b*sqrt(f)
# r = (e*d + b*sqrt(f)) / (a*sqrt(f) - e)

numerator = e * d + b * math.sqrt(f)
denominator = a * math.sqrt(f) - e

r0 = numerator / denominator

# We print the numbers that appear in the final equation that we solve:
# (3*r0 - 37) / (r0 + 4) = 1 / sqrt(2)
print("The final equation to solve for r_0 is derived from T_2(r_0) = 1/sqrt(2), which is:")
print(f"({a}*r_0 - {b}) / (r_0 + {d}) = {e}/sqrt({f})")
print("\nThe numbers in this final equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"d = {d}")
print(f"e = {e}")
print(f"f = {f}")
print("\nSolving for r_0 gives:")
print(f"r_0 = ({e}*{d} + {b}*sqrt({f})) / ({a}*sqrt({f}) - {e})")
print(f"r_0 = {r0}")

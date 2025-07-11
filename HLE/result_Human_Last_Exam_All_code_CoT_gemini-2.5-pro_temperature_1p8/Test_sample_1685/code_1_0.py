import math

# Calculate the coefficient A from the relation A^3 = -6
# We take the real cube root.
A = -math.pow(6, 1/3)

# The power n was determined to be -1 from the dominant balance analysis.
n = -1.0

# Round the coefficient to two decimal places as requested.
A_rounded = round(A, 2)
n_rounded = round(n, 2)

print("The analytical expression is of the form: y(x) = A * x^n")
print("From the analysis, the determined parameters are:")
print(f"A = {A_rounded:.2f}")
print(f"n = {n_rounded:.2f}")
print("Thus, the approximate solution for large x is:")
print(f"y(x) = {A_rounded:.2f} * x^({n_rounded:.2f})")
print("Final equation with numbers:")
print(f"y(x) = {A_rounded:.2f} / x")
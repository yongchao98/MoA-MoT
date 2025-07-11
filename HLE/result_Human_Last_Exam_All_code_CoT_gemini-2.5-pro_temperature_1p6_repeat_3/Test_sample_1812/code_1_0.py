import numpy as np

# The problem provides the following matrix A:
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# For the formula det(A) = a(ei - fh) - b(di - fg) + c(dh - eg),
# we extract the elements from the matrix.
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the intermediate values for the cofactors.
cofactor1 = (e * i - f * h)
cofactor2 = (d * i - f * g)
cofactor3 = (d * h - e * g)

# Calculate the final determinant.
determinant = a * cofactor1 - b * cofactor2 + c * cofactor3

# Print the calculation steps, showing how each number is used.
print(f"The given matrix A is:")
print(A)
print("\nThe formula for the determinant is: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

print("\nStep 1: Plugging in the numbers from matrix A:")
print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")

print("\nStep 2: Evaluating the expressions inside the parentheses:")
print(f"det(A) = {a}({e*i} - {f*h}) - ({b})({d*i} - {f*g}) + ({c})({d*h} - {e*g})")

print("\nStep 3: Simplifying the cofactors:")
print(f"det(A) = {a}({cofactor1}) - ({b})({cofactor2}) + ({c})({cofactor3})")

print("\nStep 4: Multiplying the cofactors by the row elements:")
print(f"det(A) = ({a * cofactor1}) - ({b * cofactor2}) + ({c * cofactor3})")

print("\nStep 5: Performing the final addition and subtraction:")
print(f"det(A) = {determinant}")
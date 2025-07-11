import numpy as np

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

print("The given matrix A is:")
print(A)
print("\n")

# The formula for the determinant of a 3x3 matrix [a,b,c; d,e,f; g,h,i] is:
# a(ei - fh) - b(di - fg) + c(dh - eg)
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Print the calculation steps
print("Calculating the determinant using the formula:")
print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")
term1 = a * (e * i - f * h)
term2 = b * (d * i - f * g)
term3 = c * (d * h - e * g)
print(f"det(A) = {term1} - ({term2}) + ({term3})")
det_manual = term1 - term2 + term3
print(f"det(A) = {term1} + {-term2} + {term3} = {det_manual}")
print("\n")

# Calculate the determinant using numpy
det_A = np.linalg.det(A)

print("The determinant of the matrix A is:")
# The result from numpy might be a very small floating point number close to zero
# due to floating point inaccuracies, so we can round it.
print(round(det_A))
import numpy as np

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

print("The given matrix A is:")
print(A)
print("\nTo compute the determinant of a 3x3 matrix:")
print("[[a, b, c],")
print(" [d, e, f],")
print(" [g, h, i]]")
print("We use the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")

# Extract the elements of the matrix for the formula
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Show the formula with the numbers from matrix A
print("Substituting the values from matrix A into the formula:")
calc_str = f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))"
print(calc_str)

# Calculate each term
term1_val = e*i - f*h
term2_val = d*i - f*g
term3_val = d*h - e*g
print(f"det(A) = {a}({term1_val}) - ({b})({term2_val}) + ({c})({term3_val})")

# Calculate the final result
final_result = a*term1_val - b*term2_val + c*term3_val
print(f"det(A) = {a*term1_val} - ({-b*term2_val}) + ({c*term3_val})")
print(f"det(A) = {a*term1_val} + {b*term2_val} + {c*term3_val}")
print(f"det(A) = {final_result}")

# Verify with numpy's function and print the final result
determinant = np.linalg.det(A)
print(f"\nThe final computed determinant is: {determinant}")

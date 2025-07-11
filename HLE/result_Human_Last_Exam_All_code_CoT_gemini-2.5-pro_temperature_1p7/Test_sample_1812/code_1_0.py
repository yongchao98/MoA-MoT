import numpy as np

# The problem asks to compute the determinant of a given matrix A.
# The matrix is the adjacency matrix of the Markov quiver.
# A = [[2, -1, -1],
#      [-1, 2, -1],
#      [-1, -1, 2]]
# The information m_1/2 = 13 is context and not required for the calculation.

# Define the matrix A
A = np.array([[2, -1, -1],
              [-1, 2, -1],
              [-1, -1, 2]])

# Print the matrix
print("The given adjacency matrix A is:")
print(A)
print("\nCalculating the determinant of A step by step.")
print("The formula for the determinant of a 3x3 matrix [a,b,c; d,e,f; g,h,i] is: a(e*i - f*h) - b(d*i - f*g) + c(d*h - e*g)")

# Extract elements for clarity in the printout
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Step 1: Substitute the numbers into the formula
print(f"\nStep 1: Substitute numbers from A into the formula")
print(f"det(A) = {a} * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Step 2: Perform calculations within the parentheses
print(f"\nStep 2: Calculate the products inside the parentheses")
print(f"det(A) = {a} * ({e*i} - {f*h}) - ({b}) * ({d*i} - {f*g}) + ({c}) * ({d*h} - {e*g})")

# Step 3: Simplify the expressions inside the parentheses
print(f"\nStep 3: Simplify the subtraction inside the parentheses")
print(f"det(A) = {a} * ({e*i - f*h}) - ({b}) * ({d*i - f*g}) + ({c}) * ({d*h - e*g})")

# Step 4: Perform the multiplication
print(f"\nStep 4: Perform the multiplications")
print(f"det(A) = {a * (e*i - f*h)} - {b * (d*i - f*g)} + {c * (d*h - e*g)}")

# Step 5: Final addition and subtraction
print(f"\nStep 5: Perform the final addition and subtraction")
final_det = a * (e*i - f*h) - b * (d*i - f*g) + c * (d*h - e*g)
print(f"det(A) = {a * (e*i - f*h)} - {b * (d*i - f*g)} - {abs(c * (d*h - e*g))}")
print(f"det(A) = {final_det}")


# Final Answer
print(f"\nThe determinant of the matrix is: {final_det}")

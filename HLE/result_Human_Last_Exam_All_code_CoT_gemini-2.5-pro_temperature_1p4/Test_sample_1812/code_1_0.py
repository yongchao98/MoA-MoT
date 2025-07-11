import numpy as np

# The problem asks to compute the determinant of the given adjacency matrix A.
# The matrix A is provided as:
# A = [[ 2, -1, -1],
#      [-1,  2, -1],
#      [-1, -1,  2]]
#
# The value m_{1/2} = 13 is context about the Markov quiver but not needed
# for the determinant calculation of the given matrix.
#
# We will calculate the determinant of a 3x3 matrix using Sarrus' rule:
# det(A) = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h

# Define the matrix
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract elements for clarity in the formula
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the six terms for Sarrus' rule
# Positive diagonal terms
term1 = a * e * i
term2 = b * f * g
term3 = c * d * h
# Negative diagonal terms
term4 = c * e * g
term5 = b * d * i
term6 = a * f * h

# Calculate the final determinant
determinant = term1 + term2 + term3 - term4 - term5 - term6

# Print the final equation with each number
print("The determinant is calculated using Sarrus' rule.")
print("The formula is: det(A) = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h")
print("\nSubstituting the values from the matrix A:")
print(f"det(A) = ({a})*({e})*({i}) + ({b})*({f})*({g}) + ({c})*({d})*({h}) - ({c})*({e})*({g}) - ({b})*({d})*({i}) - ({a})*({f})*({h})")
print("\nCalculating each term:")
print(f"det(A) = ({term1}) + ({term2}) + ({term3}) - ({term4}) - ({term5}) - ({term6})")
print("\nSumming the terms:")
print(f"det(A) = {determinant}")

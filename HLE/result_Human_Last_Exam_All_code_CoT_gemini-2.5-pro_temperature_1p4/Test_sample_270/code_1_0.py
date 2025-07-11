import numpy as np
import math

# Step 1: Define the matrices for the Dehn twists D_a and D_b.
A = np.array([[1, 1], [0, 1]], dtype=np.int64)
B = np.array([[1, 0], [1, 1]], dtype=np.int64)

# Step 2: Compute the matrix for the composition (D_a o D_b).
C = np.dot(A, B)

# Compute the 9th power of this matrix, which corresponds to (D_a o D_b)^9.
M = np.linalg.matrix_power(C, 9)

# Step 3: Extract the components p, s, and r from the resulting matrix M.
# M = [[p, q], [r, s]]
p = M[0, 0]
s = M[1, 1]
r = M[1, 0]

# Step 4: Calculate the numerator (the trace of M) and the denominator (the lower-left entry).
numerator = p + s
denominator = r

# The formula for the fractional Dehn twist coefficient is (p+s)/r.
# We will print the equation with the calculated numbers.
print(f"The matrix for (D_a o D_b)^9 is:\n{M}")
print("\nThe fractional Dehn twist coefficient is calculated as (trace of M) / (lower-left entry of M).")
print(f"\nThe final equation is: ({p} + {s}) / {r} = {numerator} / {denominator}")

# For completeness, we can simplify the fraction and provide the decimal value.
common_divisor = math.gcd(numerator, denominator)
simplified_num = numerator // common_divisor
simplified_den = denominator // common_divisor

print(f"The simplified fraction is: {simplified_num}/{simplified_den}")
print(f"The decimal value is: {numerator / denominator}")

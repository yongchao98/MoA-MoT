import numpy as np

# Set a random seed for reproducibility
np.random.seed(0)

# 1. Define the matrix for the birefringent medium (B).
# We'll model a quarter-wave plate with its fast axis horizontal.
# The Jones matrix for this is B.
delta = np.pi / 2  # 90-degree phase shift (quarter-wave)
B = np.array([[np.exp(-1j * delta / 2), 0],
              [0, np.exp(1j * delta / 2)]])

# 2. Define the transmission matrix for the random medium (T).
# A random medium will generally couple polarization components, so we'll
# model its Jones matrix (for a single point) as a random complex matrix.
# This matrix will not be diagonal.
T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)

# 3. Calculate the two possible orderings of the operations.
# Case 1: Light passes through B, then T. The combined operator is T @ B.
TB = T @ B
# Case 2: Light passes through T, then B. The combined operator is B @ T.
BT = B @ T

# 4. Print the results to show they are different.
# The final "equation" we are testing is TB = BT.

print("--- Demonstration of Non-Commutativity ---")
print("Matrix for Birefringent Medium (B):")
print(B)
print("\nMatrix for Random Medium (T):")
print(T)

print("\n--- Comparing the Order of Operations ---")

# Print the resulting matrix for T * B
print("\nResult of (T * B):")
print("Element (0,0):", TB[0,0])
print("Element (0,1):", TB[0,1])
print("Element (1,0):", TB[1,0])
print("Element (1,1):", TB[1,1])

# Print the resulting matrix for B * T
print("\nResult of (B * T):")
print("Element (0,0):", BT[0,0])
print("Element (0,1):", BT[0,1])
print("Element (1,0):", BT[1,0])
print("Element (1,1):", BT[1,1])

# Calculate and print the difference to make it explicit
difference = TB - BT
print("\nDifference (T*B - B*T):")
print(difference)

is_zero = np.allclose(difference, 0)
print(f"\nAre the results the same (T*B == B*T)? {is_zero}")

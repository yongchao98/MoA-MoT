import numpy as np

# Set the rotation angle for the optically active medium in radians
# Let's choose 45 degrees, which is pi/4 radians.
angle = np.pi / 4

# 1. Define the Jones Matrix 'B' for an optical rotator (circular birefringence).
# A rotator is represented by a real rotation matrix.
B = np.array([
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)]
])

# 2. Calculate the complex conjugate of B, which we call 'B_star'.
# Since B is real, B_star will be the same as B.
B_star = np.conjugate(B)

# 3. Calculate the product M = B * B_star using matrix multiplication.
# For the phase conjugation theory to hold, M should be the identity matrix.
M = B @ B_star

# 4. Define the identity matrix for comparison.
I = np.identity(2)

# --- Output the results ---
# The final equation is B * B* = M
print("The theory of reconstruction via phase conjugation requires that B * B* = I.")
print("Let's test this for an optically active medium (rotator).\n")

print("The Jones matrix B is:")
print(B)
print("\nIts complex conjugate B* is:")
print(B_star)
print("\nThe product M = B * B* is:")
print(M)
print("\nThe identity matrix I is:")
print(I)

# 5. Check if M is equal to the identity matrix.
is_identity = np.allclose(M, I)

print(f"\nIs the resulting matrix M equal to the identity matrix? {is_identity}")

if not is_identity:
    print("\nConclusion: Since B * B* is not the identity matrix, the theory does not hold for a system containing an optically active (circularly birefringent) medium.")
else:
    print("\nConclusion: B * B* is the identity matrix, so the theory holds.")

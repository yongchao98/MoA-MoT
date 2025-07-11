import numpy as np

def Rx(degrees):
    """Computes the rotation matrix around the X-axis."""
    theta = np.deg2rad(degrees)
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])

def Ry(degrees):
    """Computes the rotation matrix around the Y-axis."""
    theta = np.deg2rad(degrees)
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

def Rz(degrees):
    """Computes the rotation matrix around the Z-axis."""
    theta = np.deg2rad(degrees)
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])

# 1. Define the initial Tait-Bryan rotation parameters
alpha_tb, beta_tb, gamma_tb = 10.0, 10.0, 10.0

# Calculate the reference rotation matrix for extrinsic X-Y-Z rotation
# R_ref = Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma_tb) @ Ry(beta_tb) @ Rx(alpha_tb)

print("Reference Matrix from Tait-Bryan X(10)Y(10)Z(10) rotation:")
print(R_ref)
print("-" * 50)

# 2. Define the proper Euler angles to be tested
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# 3. Test the ZYZ convention (Answer E)
# The rotation equation for the ZYZ convention is R = Rz(alpha') * Ry(beta') * Rz(gamma')
print("Testing the ZYZ convention with the given Euler angles.")
print(f"The equation is: R_ZYZ = Rz({alpha_p}°) * Ry({beta_p}°) * Rz({gamma_p}°)")

R_ZYZ = Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p)

print("\nResulting Matrix from ZYZ convention:")
print(R_ZYZ)
print("-" * 50)

# 4. Compare the two matrices
are_matrices_equal = np.allclose(R_ref, R_ZYZ)

print(f"Are the two matrices equivalent? {are_matrices_equal}")

if are_matrices_equal:
    print("\nThe ZYZ convention (E) provides the equivalent rotation.")
else:
    print("\nThe ZYZ convention does not provide the equivalent rotation.")

<<<E>>>
import numpy as np

def Rx(degrees):
    """Creates a rotation matrix for a rotation around the X-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(degrees):
    """Creates a rotation matrix for a rotation around the Y-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(degrees):
    """Creates a rotation matrix for a rotation around the Z-axis."""
    rad = np.deg2rad(degrees)
    c, s = np.cos(rad), np.sin(rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# Step 1: Calculate the reference rotation matrix from the extrinsic Tait-Bryan angles.
alpha, beta, gamma = 10.0, 10.0, 10.0
# For extrinsic X-Y-Z, the matrix is Rz @ Ry @ Rx
R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)

print(f"The reference rotation is given by the extrinsic Tait-Bryan sequence X({alpha}°) Y({beta}°) Z({gamma}°).")
print(f"The resulting rotation matrix is R_z({gamma}) @ R_y({beta}) @ R_x({alpha}):")
print(R_ref)
print("-" * 30)

# Step 2: Define the Euler angles to test.
alpha_prime, beta_prime, gamma_prime = 139.13, 14.11, -141.05

# Step 3 & 4: Iterate through conventions, calculate their matrices, and compare.
# The conventions are given as answer choices. We only test the proper Euler conventions.
conventions = {
    "A. XZX": (Rx, Rz, Rx),
    "C. YXY": (Ry, Rx, Ry),
    "D. YZY": (Ry, Rz, Ry),
    "E. ZYZ": (Rz, Ry, Rz),
    "F. ZXZ": (Rz, Rx, Rz),
}

found_convention = None
for name, (R1, R2, R3) in conventions.items():
    # Test the intrinsic interpretation: R = R1(alpha') @ R2(beta') @ R3(gamma')
    R_candidate_intrinsic = R1(alpha_prime) @ R2(beta_prime) @ R3(gamma_prime)
    if np.allclose(R_ref, R_candidate_intrinsic, atol=1e-5):
        found_convention = name
        print(f"Found a match with the INTRINSIC convention: {name}")
        print(f"The equation is R = R1({alpha_prime}°) @ R2({beta_prime}°) @ R3({gamma_prime}°)")
        print("The calculated matrix is:")
        print(R_candidate_intrinsic)
        print("\nThis matrix matches the reference matrix.")
        break

    # Test the extrinsic interpretation: R = R3(gamma') @ R2(beta') @ R1(alpha')
    R_candidate_extrinsic = R3(gamma_prime) @ R2(beta_prime) @ R1(alpha_prime)
    if np.allclose(R_ref, R_candidate_extrinsic, atol=1e-5):
        found_convention = name
        print(f"Found a match with the EXTRINSIC convention: {name}")
        print(f"The equation is R = R3({gamma_prime}°) @ R2({beta_prime}°) @ R1({alpha_prime}°)")
        print("The calculated matrix is:")
        print(R_candidate_extrinsic)
        print("\nThis matrix matches the reference matrix.")
        break

if not found_convention:
    print("No matching convention was found among the choices.")

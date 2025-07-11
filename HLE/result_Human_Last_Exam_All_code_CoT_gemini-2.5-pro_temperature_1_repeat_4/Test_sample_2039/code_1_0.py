import numpy as np

def Rx(degrees):
    """Rotation matrix around X-axis"""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(degrees):
    """Rotation matrix around Y-axis"""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(degrees):
    """Rotation matrix around Z-axis"""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# 1. Calculate the reference rotation matrix from the extrinsic Tait-Bryan angles
alpha, beta, gamma = 10.0, 10.0, 10.0
# For extrinsic X-Y-Z, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)

print(f"The reference rotation matrix for the extrinsic X({alpha}°) Y({beta}°) Z({gamma}°) rotation is:")
print(R_ref)
print("-" * 30)

# 2. Define the equivalent Euler angles
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# 3. Test each proper Euler convention (assuming intrinsic rotations)
conventions = {
    "A. XZX": Rx(alpha_p) @ Rz(beta_p) @ Rx(gamma_p),
    "C. YXY": Ry(alpha_p) @ Rx(beta_p) @ Ry(gamma_p),
    "D. YZY": Ry(alpha_p) @ Rz(beta_p) @ Ry(gamma_p),
    "E. ZYZ": Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p),
    "F. ZXZ": Rz(alpha_p) @ Rx(beta_p) @ Rz(gamma_p),
}

# 4. Find and print the matching convention and its matrix
print(f"Testing Euler angles a'={alpha_p}°, b'={beta_p}°, g'={gamma_p}° for each convention:")
for name, R_candidate in conventions.items():
    # Using a small tolerance to account for floating point inaccuracies
    if np.allclose(R_ref, R_candidate, atol=1e-5):
        print(f"\nFound a match with convention: {name}")
        print("The final equation is the equality between the two matrices:")
        print(f"R_ref = R_candidate for convention {name}")
        
        print("\nMatrix from initial rotation (R_ref):")
        # Outputting each number in the matrix
        for row in R_ref:
            print(f"[{row[0]:.8f} {row[1]:.8f} {row[2]:.8f}]")
            
        print(f"\nMatrix from {name} convention (R_candidate):")
        # Outputting each number in the matrix
        for row in R_candidate:
            print(f"[{row[0]:.8f} {row[1]:.8f} {row[2]:.8f}]")
        
        print(f"\nThe correct convention is {name.split('.')[0]}.")
        break
<<<F>>>
import numpy as np

def Rx(theta_deg):
    """Rotation matrix around X-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta_deg):
    """Rotation matrix around Y-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(theta_deg):
    """Rotation matrix around Z-axis"""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# 1. Define initial extrinsic rotation angles.
alpha, beta, gamma = 10.0, 10.0, 10.0

# 2. Calculate the reference rotation matrix from the extrinsic XYZ rotation.
# The order of matrix multiplication is Rz * Ry * Rx for extrinsic XYZ.
R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# 3. Define the equivalent proper Euler angles.
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# 4. Define the possible proper Euler angle conventions and their matrix functions.
# The conventions are for intrinsic rotations, so the matrix order is sequential.
conventions = {
    "A. XZX": (Rx, Rz, Rx),
    "C. YXY": (Ry, Rx, Ry),
    "D. YZY": (Ry, Rz, Ry),
    "E. ZYZ": (Rz, Ry, Rz),
    "F. ZXZ": (Rz, Rx, Rz),
}

found_convention = None
# 5. Test each convention.
for name, (R1, R2, R3) in conventions.items():
    R_test = R1(alpha_p) @ R2(beta_p) @ R3(gamma_p)
    if np.allclose(R_ref, R_test, atol=1e-5):
        found_convention = name
        
        # 6. Print the details for the matching convention.
        print(f"The initial extrinsic rotation is given by X({alpha} deg), Y({beta} deg), Z({gamma} deg).")
        print(f"The reference rotation matrix R_ref = R_z({gamma}) @ R_y({beta}) @ R_x({alpha}) is:")
        print(R_ref)
        print("\n" + "="*50 + "\n")
        print(f"Testing the '{name}' convention with angles alpha'={alpha_p}, beta'={beta_p}, gamma'={gamma_p}.")
        convention_formula = f"R_test = R_{name[3]}({alpha_p}) @ R_{name[4]}({beta_p}) @ R_{name[5]}({gamma_p})"
        print(f"The formula is: {convention_formula}")
        print("The resulting test matrix R_test is:")
        print(R_test)
        print("\n" + "="*50 + "\n")
        print("Comparing R_ref and R_test shows they are equivalent.")
        print(f"Therefore, the correct convention is {name}.")
        break

if not found_convention:
    print("No matching convention was found among the proper Euler angle choices.")

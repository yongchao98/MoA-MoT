import numpy as np

def Rx(theta_deg):
    """Creates a rotation matrix for a rotation around the X-axis."""
    theta_rad = np.radians(theta_deg)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(theta_deg):
    """Creates a rotation matrix for a rotation around the Y-axis."""
    theta_rad = np.radians(theta_deg)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(theta_deg):
    """Creates a rotation matrix for a rotation around the Z-axis."""
    theta_rad = np.radians(theta_deg)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# Step 1: Define the angles for the initial extrinsic rotation
alpha, beta, gamma = 10.0, 10.0, 10.0

# Step 2: Calculate the target rotation matrix from the extrinsic X-Y-Z rotation.
# The order of multiplication is Rz * Ry * Rx.
R_target = Rz(gamma) @ Ry(beta) @ Rx(alpha)

print("The target rotation matrix from the extrinsic X(10)Y(10)Z(10) rotation is:")
print(R_target)
print("-" * 60)

# Step 3: Define the given Euler angles to be tested
alpha_prime, beta_prime, gamma_prime = 139.13, 14.11, -141.05

# Step 4: Define the conventions to be tested
# For intrinsic rotations (like Euler angles), the matrix order is R_A * R_B * R_C.
conventions_to_test = {
    "A. XZX": Rx(alpha_prime) @ Rz(beta_prime) @ Rx(gamma_prime),
    "B. XYZ": Rx(alpha_prime) @ Ry(beta_prime) @ Rz(gamma_prime),
    "C. YXY": Ry(alpha_prime) @ Rx(beta_prime) @ Ry(gamma_prime),
    "D. YZY": Ry(alpha_prime) @ Rz(beta_prime) @ Ry(gamma_prime),
    "E. ZYZ": Rz(alpha_prime) @ Ry(beta_prime) @ Rz(gamma_prime),
    "F. ZXZ": Rz(alpha_prime) @ Rx(beta_prime) @ Rz(gamma_prime)
}

# Step 5 & 6: Iterate through conventions, find the match, and print the result
found = False
for name, R_test in conventions_to_test.items():
    if np.allclose(R_target, R_test):
        print(f"Found a match with convention: {name}")
        print("\nThe final equation demonstrates the two matrices are equivalent:")
        print(f"\nR_target = Rz({gamma}) @ Ry({beta}) @ Rx({alpha}) = \n{R_target}\n")
        
        convention_name = name.split(". ")[1]
        alpha_str = f"R{convention_name[0]}({alpha_prime})"
        beta_str = f"R{convention_name[1]}({beta_prime})"
        gamma_str = f"R{convention_name[2]}({gamma_prime})"
        
        print(f"R_test = {alpha_str} @ {beta_str} @ {gamma_str} = \n{R_test}\n")
        
        print(f"The correct convention is {name.split('.')[0]}.")
        found = True
        break

if not found:
    print("No matching convention was found among the options.")

<<<F>>>
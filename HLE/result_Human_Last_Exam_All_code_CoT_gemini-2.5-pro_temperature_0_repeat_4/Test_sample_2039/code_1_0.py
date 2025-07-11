import numpy as np

def Rx(theta_deg):
    """Rotation matrix around the X-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta_deg):
    """Rotation matrix around the Y-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(theta_deg):
    """Rotation matrix around the Z-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# 1. Define the initial Tait-Bryan rotation angles
alpha = 10.0
beta = 10.0
gamma = 10.0

# Calculate the initial rotation matrix for extrinsic X-Y-Z convention
# R_initial = Rz(gamma) @ Ry(beta) @ Rx(alpha)
R_initial = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# 2. Define the equivalent proper Euler angles
alpha_prime = 139.13
beta_prime = 14.11
gamma_prime = -141.05

# 3. Define the possible conventions and their corresponding matrix multiplications
# The key is the answer choice, the value is a tuple with (name, R_axis1, R_axis2, R_axis3)
conventions = {
    "A": ("XZX", Rx, Rz, Rx),
    "B": ("XYZ", Rx, Ry, Rz),
    "C": ("YXY", Ry, Rx, Ry),
    "D": ("YZY", Ry, Rz, Ry),
    "E": ("ZYZ", Rz, Ry, Rz),
    "F": ("ZXZ", Rz, Rx, Rz),
}

# 4. Iterate through conventions to find the match
print("Searching for the correct Euler angle convention...")
found_match = False
for key, (name, R1, R2, R3) in conventions.items():
    # Calculate the test matrix for the current convention (intrinsic rotation)
    R_test = R1(alpha_prime) @ R2(beta_prime) @ R3(gamma_prime)

    # Check if the test matrix is numerically close to the initial matrix
    if np.allclose(R_initial, R_test):
        print(f"\nMatch found! The correct convention is {key}: {name}")
        
        print("\nThe initial rotation from extrinsic X-Y-Z is equivalent to the rotation from the found convention.")
        print("The equality of the rotation matrices is shown below:")
        
        print(f"\nR_Z({gamma}) @ R_Y({beta}) @ R_X({alpha}) =")
        print(R_initial)
        
        print(f"\n{name[0]}_R({alpha_prime}) @ {name[1]}_R({beta_prime}) @ {name[2]}_R({gamma_prime}) =")
        print(R_test)
        
        print("\nTherefore, the final equation with all the numbers is:")
        print(f"R_Z({gamma}) R_Y({beta}) R_X({alpha}) = R_Z({alpha_prime}) R_X({beta_prime}) R_Z({gamma_prime})")
        
        found_match = True
        break

if not found_match:
    print("\nNo matching convention was found among the options.")

<<<F>>>
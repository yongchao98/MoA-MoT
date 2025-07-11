import numpy as np

def Rx(theta_deg):
    """Rotation matrix around X-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta_rad), -np.sin(theta_rad)],
        [0, np.sin(theta_rad), np.cos(theta_rad)]
    ])

def Ry(theta_deg):
    """Rotation matrix around Y-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    return np.array([
        [np.cos(theta_rad), 0, np.sin(theta_rad)],
        [0, 1, 0],
        [-np.sin(theta_rad), 0, np.cos(theta_rad)]
    ])

def Rz(theta_deg):
    """Rotation matrix around Z-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    return np.array([
        [np.cos(theta_rad), -np.sin(theta_rad), 0],
        [np.sin(theta_rad), np.cos(theta_rad), 0],
        [0, 0, 1]
    ])

# 1. Calculate the reference rotation matrix from Tait-Bryan angles
alpha, beta, gamma = 10.0, 10.0, 10.0
# The convention is extrinsic X-Y-Z, so R_ref = Rx(alpha) @ Ry(beta) @ Rz(gamma)
R_ref = Rx(alpha) @ Ry(beta) @ Rz(gamma)

# 2. Define the Euler angles to test
ap, bp, cp = 139.13, 14.11, -141.05

# 3. Define the conventions to test
conventions = {
    'A': 'XZX',
    'C': 'YXY',
    'D': 'YZY',
    'E': 'ZYZ',
    'F': 'ZXZ'
}
rot_funcs = {'X': Rx, 'Y': Ry, 'Z': Rz}

# 4. Iterate through conventions and find the match
found = False
for key, conv_str in conventions.items():
    axis1_char, axis2_char, axis3_char = conv_str[0], conv_str[1], conv_str[2]
    
    R1 = rot_funcs[axis1_char]
    R2 = rot_funcs[axis2_char]
    R3 = rot_funcs[axis3_char]
    
    # Calculate the test matrix R_test = R_axis1(ap) @ R_axis2(bp) @ R_axis3(cp)
    R_test = R1(ap) @ R2(bp) @ R3(cp)
    
    # Compare with reference matrix using a tolerance for floating point numbers
    if np.allclose(R_ref, R_test, atol=1e-4):
        print(f"The correct Euler angle convention is: {conv_str} (Choice {key})")
        print("This convention establishes the following equivalence:")
        print(f"R_X({alpha}) @ R_Y({beta}) @ R_Z({gamma}) is equivalent to R_{axis1_char}({ap}) @ R_{axis2_char}({bp}) @ R_{axis3_char}({cp})")
        found = True
        break

if not found:
    print("No matching convention was found among the given options.")

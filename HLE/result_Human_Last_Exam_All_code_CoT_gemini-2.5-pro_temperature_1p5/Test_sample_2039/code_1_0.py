import numpy as np

def Rx(theta_deg):
    """Creates a 3x3 rotation matrix around the X-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(theta_deg):
    """Creates a 3x3 rotation matrix around the Y-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(theta_deg):
    """Creates a 3x3 rotation matrix around the Z-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# 1. Define Tait-Bryan angles and calculate the target rotation matrix.
# The convention is extrinsic X-Y-Z, so the matrix is Rz*Ry*Rx.
alpha, beta, gamma = 10.0, 10.0, 10.0
R_target = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# 2. Define the equivalent Euler angles to be tested.
ap, bp, gp = 139.13, 14.11, -141.05

# 3. Define the possible proper Euler conventions (intrinsic rotations).
# Each entry is a tuple: (Answer Choice, Convention Name, [Rotation Functions])
conventions_to_test = [
    ('A', 'XZX', [Rx, Rz, Rx]),
    ('C', 'YXY', [Ry, Rx, Ry]),
    ('D', 'YZY', [Ry, Rz, Ry]),
    ('E', 'ZYZ', [Rz, Ry, Rz]),
    ('F', 'ZXZ', [Rz, Rx, Rz]),
]

# 4. Iterate through conventions, calculate candidate matrix, and check for a match.
found_answer = None
for key, name, funcs in conventions_to_test:
    R1, R2, R3 = funcs
    R_candidate = R1(ap) @ R2(bp) @ R3(gp)
    
    # Use np.allclose for robust floating-point comparison
    if np.allclose(R_target, R_candidate, atol=1e-4):
        found_answer = key
        print(f"Found a match with the {name} convention.")
        print("The equivalent rotation can be expressed as:")
        # Print the final equation with all the numbers
        print(f"R_extrinsic-XYZ({alpha:.2f}°, {beta:.2f}°, {gamma:.2f}°) = R_intrinsic-{name}({ap:.2f}°, {bp:.2f}°, {gp:.2f}°)")
        break

if not found_answer:
    print("No matching convention was found among the options.")

# Final answer in the specified format
if found_answer:
    print(f"\n<<<__{found_answer}__>>>")

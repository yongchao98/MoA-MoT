import numpy as np

def Rx(theta_deg):
    """Rotation matrix around the X-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(theta_deg):
    """Rotation matrix around the Y-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(theta_deg):
    """Rotation matrix around the Z-axis."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# 1. Define initial Tait-Bryan angles and calculate the reference matrix
tb_alpha, tb_beta, tb_gamma = 10, 10, 10
# For extrinsic X-Y-Z rotation, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(tb_gamma) @ Ry(tb_beta) @ Rx(tb_alpha)

# 2. Define the equivalent proper Euler angles
e_alpha, e_beta, e_gamma = 139.13, 14.11, -141.05

# 3. Define the conventions to test from the answer choices
# For an extrinsic convention 'ABC' with angles a,b,c the matrix is Rc(c) @ Rb(b) @ Ra(a)
conventions = {
    "A": {"name": "XZX", "matrix": Rx(e_gamma) @ Rz(e_beta) @ Rx(e_alpha)},
    "C": {"name": "YXY", "matrix": Ry(e_gamma) @ Rx(e_beta) @ Ry(e_alpha)},
    "D": {"name": "YZY", "matrix": Ry(e_gamma) @ Rz(e_beta) @ Ry(e_alpha)},
    "E": {"name": "ZYZ", "matrix": Rz(e_gamma) @ Ry(e_beta) @ Rz(e_alpha)},
    "F": {"name": "ZXZ", "matrix": Rz(e_gamma) @ Rx(e_beta) @ Rz(e_alpha)},
}

# 4. Find the matching convention
found_answer = None
for key, value in conventions.items():
    if np.allclose(R_ref, value["matrix"], atol=1e-4):
        found_answer = key
        print(f"Match Found: Convention {value['name']} ({key})")
        print("The final rotation is given by the matrix equation:")
        print(f"R = R_{value['name'][2]}({e_gamma}°) * R_{value['name'][1]}({e_beta}°) * R_{value['name'][0]}({e_alpha}°)\n")

        print("Reference Tait-Bryan Matrix R_Z(10)R_Y(10)R_X(10):")
        print(R_ref)
        print(f"\nEquivalent {value['name']} Matrix:")
        print(value["matrix"])
        break

if found_answer is None:
    print("No matching convention found among the choices.")

print(f"\n<<<E>>>")
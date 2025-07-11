import numpy as np

def Rx(theta_deg):
    """Rotation matrix around X-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta_deg):
    """Rotation matrix around Y-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def Rz(theta_deg):
    """Rotation matrix around Z-axis for an angle in degrees."""
    theta_rad = np.deg2rad(theta_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

# --- Step 1: Calculate the original rotation matrix ---
# Tait-Bryan angles
alpha, beta, gamma = 10.0, 10.0, 10.0
print(f"Calculating reference rotation matrix for extrinsic X-Y-Z rotation with:")
print(f"alpha = {alpha} deg, beta = {beta} deg, gamma = {gamma} deg\n")

# Extrinsic X-Y-Z rotation: R = Rz(gamma) * Ry(beta) * Rx(alpha)
R_orig = Rz(gamma) @ Ry(beta) @ Rx(alpha)

# --- Step 2: Define Euler angles and conventions to test ---
# Equivalent proper Euler angles
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05
print(f"Testing proper Euler angle conventions with:")
print(f"alpha' = {alpha_p} deg, beta' = {beta_p} deg, gamma' = {gamma_p} deg\n")


# Conventions to test (assuming intrinsic rotations)
# For an intrinsic convention ABC with angles a,b,c, R = Ra(a) * Rb(b) * Rc(c)
conventions = {
    "A": ("XZX", Rx(alpha_p) @ Rz(beta_p) @ Rx(gamma_p)),
    "C": ("YXY", Ry(alpha_p) @ Rx(beta_p) @ Ry(gamma_p)),
    "D": ("YZY", Ry(alpha_p) @ Rz(beta_p) @ Ry(gamma_p)),
    "E": ("ZYZ", Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p)),
    "F": ("ZXZ", Rz(alpha_p) @ Rx(beta_p) @ Rz(gamma_p)),
}

# --- Step 3: Find the matching convention ---
found_match = False
for key, (name, R_test) in conventions.items():
    # Use a tolerance to account for rounding in the problem's given angles
    if np.allclose(R_orig, R_test, atol=1e-3):
        print(f"Found a match: Convention {key} ({name})")
        print("Original Matrix (from Tait-Bryan angles):")
        print(R_orig)
        print(f"\nTest Matrix (from {name} Euler angles):")
        print(R_test)
        print("\nThe matrices are equivalent.")
        final_answer = key
        found_match = True
        break

if not found_match:
    print("No matching convention was found.")

print(f"\n<<<F>>>")